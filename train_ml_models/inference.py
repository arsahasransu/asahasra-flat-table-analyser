import awkward as ak
import numpy as np
import torch
import uproot

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from ml_helper_functions import SoftIsoSumNetwork, tkel_bs, puppi_bs, puppicands
from ml_helper_functions import split_puppi_by_pdgid, flatten_puppi_collections
from ml_helper_functions import flatten_tkel_collection, pair_tkel_puppi_and_filter
from ml_helper_functions import get_tkel_branches, get_puppi_branches

def add_soft_iso_to_root_file(file_path, tkeltag, puppitag, *,
                              tree_name='snapshot', 
                              model_dir='train_ml_models', 
                              model_suffix='', batch_size=1024):
    """
    Reads a ROOT file using uproot, computes the Soft IsoSum using the trained PyTorch model,
    and updates the ROOT file by rewriting the tree with the new 'soft_iso_sum' branch added.
    
    Parameters:
      file_path: path to the *_snapshot.root file.
      tkeltag: tag for TkEle branches (e.g., 'Pt5_EB_MCH_')
      puppitag: tag for Puppi branches (e.g., 'Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_')
      tree_name: name of the tree in the ROOT file (default: 'snapshot')
      model_dir: path to the directory containing `best_soft_iso_model.pt` and `norm_stats.pt`.
      model_suffix: suffix for model and stat filenames, e.g., '_eb' or '_ee'.
    """
    tkel_branches = get_tkel_branches(tkeltag)
    puppi_branches = get_puppi_branches(puppitag)
    
    with uproot.open(file_path) as f:
        tree = f[tree_name]
        data = tree.arrays()
        
    # Build tkel and puppi records
    tkel = ak.zip({k: data[v] for k, v in tkel_branches.items()})
    puppi = ak.zip({k: data[v] for k, v in puppi_branches.items()})
    
    # 2. Cartesian product per event: N_events x N_tkels x N_puppi
    tkel_with_puppi = pair_tkel_puppi_and_filter(tkel, puppi, dr_cut=0.5)
    
    # 3. Split by pdgId
    tkel_with_puppi = split_puppi_by_pdgid(tkel_with_puppi)
    
    # 4. Flatten puppi collections
    tkel_with_puppi = flatten_puppi_collections(tkel_with_puppi, max_items=5)
    
    # 5. Flatten tkel collection
    flat_recs = flatten_tkel_collection(tkel_with_puppi)
    
    # Keep track of the event structure
    num_tkels = ak.num(flat_recs, axis=1)
    
    # Flatten across events to get a 1D array of tkels
    per_tkel_data = ak.flatten(flat_recs, axis=1)
    
    # Ensure correct order of features to match what was trained
    # We load the training features to get the exact order
    
    train_data_path = os.path.join(model_dir, 'convert_root_to_torch', f'{model_suffix.strip("_")}_train_data.pt')
    if os.path.exists(train_data_path):
        train_info = torch.load(train_data_path, weights_only=True)
        feature_names = [f for f in train_info['features'] if f != 'label']
    else:
        # Fallback to reconstructing the order from code logic
        collections = ['mpuppiel', 'mpuppiph', 'mpuppimu', 'mpuppich', 'mpuppinh']
        max_items = 5
        tkel_fields = [f"tkel_{f}" for f in ['pt', 'eta', 'phi', 'caloEta', 'caloPhi', 'tkPt', 'tkEta', 'tkPhi', 'charge', 'hwQual', 'vz']]
        puppi_fields = []
        for coll in collections:
            for i in range(max_items):
                puppi_fields.extend([f"{coll}_{i}_{f}" for f in ['pt', 'eta', 'phi', 'mass', 'charge', 'dxy', 'hwDxy', 'hwTkQuality', 'pdgId', 'puppiWeight', 'z0']])
        feature_names = tkel_fields + puppi_fields
        
    # Stack features
    features = []
    for f in feature_names:
        features.append(ak.to_numpy(per_tkel_data[f]))
    
    X = np.stack(features, axis=-1)
    X_tensor = torch.tensor(X, dtype=torch.float32)
    
    tkel_tensor = X_tensor[:, 0:tkel_bs]
    puppi_tensor_flat = X_tensor[:, tkel_bs:]
    puppi_tensor = puppi_tensor_flat.reshape(-1, puppicands, puppi_bs)
    puppi_pt = puppi_tensor[:, :, 0].clone()
    
    # Apply log(pT) transform
    tkel_tensor[:, 0] = torch.log(torch.clamp(tkel_tensor[:, 0], min=1e-6))
    puppi_tensor[:, :, 0] = torch.log(torch.clamp(puppi_tensor[:, :, 0], min=1e-6))
    
    # Apply normalization
    norm_stats = torch.load(os.path.join(model_dir, f'norm_stats{model_suffix}.pt'), weights_only=True)
    tkel_norm = (tkel_tensor - norm_stats['tkel_mean']) / norm_stats['tkel_std']
    puppi_norm = (puppi_tensor - norm_stats['puppi_mean']) / norm_stats['puppi_std']
    
    # Run inference
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = SoftIsoSumNetwork(hidden_dims=[32, 32]).to(device)
    model.load_state_dict(torch.load(os.path.join(model_dir, f'best_soft_iso_model{model_suffix}.pt'), map_location=device, weights_only=True))
    model.eval()
    
    iso_sums = []
    dataset = torch.utils.data.TensorDataset(tkel_norm, puppi_norm, puppi_pt)
    loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False)
    
    with torch.no_grad():
        for t_n, p_n, p_pt in loader:
            t_n, p_n, p_pt = t_n.to(device), p_n.to(device), p_pt.to(device)
            _, batch_iso_sum, _ = model(t_n, p_n, p_pt)
            iso_sums.append(batch_iso_sum.cpu().numpy())
            
    iso_sum_flat = np.concatenate(iso_sums)
    
    # Unflatten back to event structure (N_events x N_tkels)
    iso_sum_unflat = ak.unflatten(iso_sum_flat, num_tkels)
    
    # Add back to data and recreate the file
    new_branch_name = f"TkEleL2_{tkeltag}weighted_iso_score"
    data[new_branch_name] = iso_sum_unflat
    
    with uproot.recreate(file_path) as f:
        f[tree_name] = data
        
    print(f"Successfully added {new_branch_name} to {file_path}")

if __name__ == '__main__':
    # Add main block to process all snapshot files in the current directory
    snapshot_files = [
        ('DY_PU200_EB_snapshot.root', 'Pt5_EB_MCH_', 'Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_', '_eb'),
        ('MinBias_EB_snapshot.root', 'Pt5_EB_', 'Pt1_TkEleL2Pt5EB_0p0dR0p5_', '_eb'),
        ('DY_PU200_EE_snapshot.root', 'Pt5_EE_MCH_', 'Pt1_TkEleL2Pt5EEMCH_0p0dR0p5_', '_ee'),
        ('MinBias_EE_snapshot.root', 'Pt5_EE_', 'Pt1_TkEleL2Pt5EE_0p0dR0p5_', '_ee')
    ]
    
    model_dir = os.path.dirname(__file__)
    
    for fname, tkeltag, puppitag, suffix in snapshot_files:
        if os.path.exists(fname):
            print(f"Processing {fname}...")
            add_soft_iso_to_root_file(fname, tkeltag, puppitag, model_dir=model_dir, model_suffix=suffix)
        else:
            print(f"File {fname} not found, skipping.")
