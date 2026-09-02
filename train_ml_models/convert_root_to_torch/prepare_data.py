import awkward as ak
import numpy as np
import uproot
import torch


import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from ml_helper_functions import get_tkel_branches, get_puppi_branches
from ml_helper_functions import split_puppi_by_pdgid, flatten_puppi_collections
from ml_helper_functions import flatten_tkel_collection, pair_tkel_puppi_and_filter

def flatten_per_tkel(filename: str, *,
                     tkeltag = "",
                     puppitag = ""):

    try:
        file = uproot.open(filename)['snapshot']
    except:
        raise FileExistsError(f'Unable to snapshot tree in file: {filename}')
    print(f'Opened {filename}')

    # List the exact branches we need to save memory and time
    tkel_branches = get_tkel_branches(tkeltag)
    puppi_branches = get_puppi_branches(puppitag)
    
    branches_to_load = list(tkel_branches.values()) + list(puppi_branches.values())
    data = file.arrays(filter_name=branches_to_load)
    print('Converted data to awkward arrays')
    
    tkel = ak.zip({k: data[v] for k, v in tkel_branches.items()})
    puppi = ak.zip({k: data[v] for k, v in puppi_branches.items()})

    tkel_with_puppi = pair_tkel_puppi_and_filter(tkel, puppi, dr_cut=0.5)

    per_tkel_data = ak.flatten(tkel_with_puppi, axis=1)
    return per_tkel_data


def save_as_pytorch(data, filename):
    feature_names = [f for f in data.fields if f not in ('label', 'weight')]
    features = []
    for f in feature_names:
        features.append(ak.to_numpy(data[f]))

    X = np.stack(features, axis=-1)
    X_tensor = torch.tensor(X, dtype=torch.float32)

    save_dict = {"x": X_tensor, "features": feature_names}

    if 'label' in data.fields:
        y_tensor = torch.tensor(ak.to_numpy(data['label']), dtype=torch.long)
        save_dict["y"] = y_tensor

    if 'weight' in data.fields:
        w_tensor = torch.tensor(ak.to_numpy(data['weight']), dtype=torch.float32)
        save_dict["w"] = w_tensor

    torch.save(save_dict, filename)
    print(f"Saved to {filename}")

def prepare_ml_data(signal_file: str, bkg_file: str, *,
                    stkeltag = "", btkeltag = "",
                    spuppitag = "", bpuppitag = "",
                    max_items = 2):
    print("Processing signal...")
    sig_data = flatten_per_tkel(signal_file,
                                tkeltag = stkeltag,
                                puppitag = spuppitag)
    if sig_data is None: return None, None
    sig_data = split_puppi_by_pdgid(sig_data)
    sig_data = flatten_puppi_collections(sig_data, max_items=max_items)
    sig_data = flatten_tkel_collection(sig_data)
    sig_data["label"] = 1
    n_sig = len(sig_data)
    print(f"Signal events: {n_sig}")
    
    print("Processing background...")
    bkg_data = flatten_per_tkel(bkg_file,
                                tkeltag = btkeltag,
                                puppitag = bpuppitag)
    if bkg_data is None: return None, None
    bkg_data = split_puppi_by_pdgid(bkg_data)
    bkg_data = flatten_puppi_collections(bkg_data, max_items=max_items)
    bkg_data = flatten_tkel_collection(bkg_data)
    bkg_data["label"] = 0
    n_bkg = len(bkg_data)
    print(f"Background events: {n_bkg}")

    # Down-weight background so that signal and background contribute equally.
    # Background weight = 1.0, signal weight = N_bkg / N_sig
    bkg_data["weight"] = ak.full_like(bkg_data["label"], 1.0)
    sig_data["weight"] = ak.full_like(sig_data["label"], float(n_bkg) / float(n_sig))
    print(f"Signal weight: {n_bkg / n_sig:.2f} (bkg/sig ratio = {n_bkg / n_sig:.2f})")
    
    print("Combining and shuffling datasets...")
    combined = ak.concatenate([sig_data, bkg_data])
    
    # Shuffle
    indices = np.random.permutation(len(combined))
    shuffled = combined[indices]
    
    # Split 80:20
    split_idx = int(0.8 * len(shuffled))
    train_data = shuffled[:split_idx]
    test_data = shuffled[split_idx:]
    
    print(f"Dataset split: {len(train_data)} training events, {len(test_data)} testing events")
    return train_data, test_data


if __name__ == "__main__":
    eb_res = prepare_ml_data(
        '../../DoubleElectronGun_PU200_EB_snapshot.root',
        '../../MinBias_EB_snapshot.root',
        stkeltag = "Pt5_EB_MCH_", btkeltag = "Pt5_EB_",
        spuppitag = "Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_",
        bpuppitag = "Pt1_TkEleL2Pt5EB_0p0dR0p5_",
        max_items = 5
    )
    if eb_res[0] is not None:
        eb_train_data, eb_test_data = eb_res
        save_as_pytorch(eb_train_data, 'eb_train_data.pt')
        save_as_pytorch(eb_test_data, 'eb_test_data.pt')

    ee_res = prepare_ml_data(
        '../../DoubleElectronGun_PU200_EE_snapshot.root',
        '../../MinBias_EE_snapshot.root',
        stkeltag = "Pt5_EE_MCH_", btkeltag = "Pt5_EE_",
        spuppitag = "Pt1_TkEleL2Pt5EEMCH_0p0dR0p5_",
        bpuppitag = "Pt1_TkEleL2Pt5EE_0p0dR0p5_",
        max_items = 5
    )
    if ee_res[0] is not None:
        ee_train_data, ee_test_data = ee_res
        save_as_pytorch(ee_train_data, 'ee_train_data.pt')
        save_as_pytorch(ee_test_data, 'ee_test_data.pt')
