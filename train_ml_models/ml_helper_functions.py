import awkward as ak
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

tkel_bs = 11
puppi_bs = 11 + 3  # dR, deta, dphi
nlargestpt_per_puppipart = 5
npuppipart = 5
puppicands = nlargestpt_per_puppipart * npuppipart

class SoftIsoSumNetwork(nn.Module):
    def __init__(self, tkel_dim=tkel_bs, puppi_dim=puppi_bs, hidden_dims=[64, 32]):
        super().__init__()
        
        # Build shared MLP for candidates
        layers = []
        input_dim = tkel_dim + puppi_dim
        for h_dim in hidden_dims:
            layers.append(nn.Linear(input_dim, h_dim))
            layers.append(nn.ReLU())
            input_dim = h_dim
        layers.append(nn.Linear(input_dim, 1))
        
        self.mlp = nn.Sequential(*layers)
        
        # Learnable threshold and scale for BCE conversion
        # We want: logit = (threshold - iso_sum) * scale
        self.threshold = nn.Parameter(torch.tensor([5.0]))
        self.scale_raw = nn.Parameter(torch.tensor([0.0])) # F.softplus(0) ~ 0.69
        
    def forward(self, tkel_norm, puppi_norm, puppi_pt_unnorm):
        """
        tkel_norm: (batch, tkel_bs)
        puppi_norm: (batch, puppicands, puppi_bs)
        puppi_pt_unnorm: (batch, puppicands) - physical pT to sum over
        """
        batch_size = tkel_norm.shape[0]
        
        # Expand tkel: (batch, puppicands, tkel_bs)
        tkel_expanded = tkel_norm.unsqueeze(1).expand(-1, puppicands, -1)
        
        # Concat: (batch, puppicands, tkel_bs + puppi_bs)
        combined = torch.cat([tkel_expanded, puppi_norm], dim=2)
        
        # MLP -> Weights: (batch, puppicands)
        weights_logit = self.mlp(combined).squeeze(-1)
        weights = torch.sigmoid(weights_logit)
        
        # Weighted sum of physical pT
        iso_sum = torch.sum(weights * puppi_pt_unnorm, dim=1) # (batch,)
        
        # Convert to BCE logit
        scale = F.softplus(self.scale_raw)
        logits = (self.threshold - iso_sum) * scale
        
        return logits, iso_sum, weights

def get_tkel_branches(tkeltag):
    return {
        "pt": f"TkEleL2_{tkeltag}pt",
        "eta": f"TkEleL2_{tkeltag}eta",
        "phi": f"TkEleL2_{tkeltag}phi",
        "caloEta": f"TkEleL2_{tkeltag}caloEta",
        "caloPhi": f"TkEleL2_{tkeltag}caloPhi",
        "tkPt": f"TkEleL2_{tkeltag}tkPt",
        "tkEta": f"TkEleL2_{tkeltag}tkEta",
        "tkPhi": f"TkEleL2_{tkeltag}tkPhi",
        "charge": f"TkEleL2_{tkeltag}charge",
        "hwQual": f"TkEleL2_{tkeltag}hwQual",
        "vz": f"TkEleL2_{tkeltag}vz"
    }

def get_puppi_branches(puppitag):
    return {
        "pt": f"L1PuppiCands_{puppitag}pt",
        "eta": f"L1PuppiCands_{puppitag}eta",
        "phi": f"L1PuppiCands_{puppitag}phi",
        "mass": f"L1PuppiCands_{puppitag}mass",
        "charge": f"L1PuppiCands_{puppitag}charge",
        "dxy": f"L1PuppiCands_{puppitag}dxy",
        "hwDxy": f"L1PuppiCands_{puppitag}hwDxy",
        "hwTkQuality": f"L1PuppiCands_{puppitag}hwTkQuality",
        "pdgId": f"L1PuppiCands_{puppitag}pdgId",
        "puppiWeight": f"L1PuppiCands_{puppitag}puppiWeight",
        "z0": f"L1PuppiCands_{puppitag}z0"
    }

def split_puppi_by_pdgid(data):
    """
    Splits the nested 'mpuppi' collection into separate collections based on absolute pdgId:
    - mpuppiel: 11
    - mpuppiph: 22
    - mpuppimu: 13
    - mpuppich: 211
    - mpuppinh: 130
    """
    mpuppi = data.mpuppi
    abs_pdg = np.abs(mpuppi.pdgId)
    
    data["mpuppiel"] = mpuppi[abs_pdg == 11]
    data["mpuppiph"] = mpuppi[abs_pdg == 22]
    data["mpuppimu"] = mpuppi[abs_pdg == 13]
    data["mpuppich"] = mpuppi[abs_pdg == 211]
    data["mpuppinh"] = mpuppi[abs_pdg == 130]
    
    return data

def flatten_puppi_collections(data, max_items=5):
    """
    Takes the max_items highest pT elements from each mpuppi category,
    pads with 0 if fewer than max_items exist, and flattens them into individual columns.
    Original mpuppi collections are removed for efficiency.
    """
    collections = ['mpuppiel', 'mpuppiph', 'mpuppimu', 'mpuppich', 'mpuppinh']
    
    new_cols = {}
    
    # Retain fields that are not part of the PUPPI collections and not the general 'mpuppi'
    for field in data.fields:
        if field not in collections and field != "mpuppi":
            new_cols[field] = data[field]
        
    for coll in collections:
        if coll not in data.fields:
            print(f"Warning, did not find {coll} in data")
            continue
            
        # Handle case where the collection is nested one level deeper (like in inference per event vs per tkel)
        # We need to sort by pt
        pt = data[coll].pt
        # If it's a 3D array (events x tkels x puppi), argsort needs axis=2
        # If it's a 2D array (tkel x puppi), argsort needs axis=1
        axis = 2 if pt.ndim == 3 else 1
        sorted_coll = data[coll][ak.argsort(pt, ascending=False, axis=axis)]
        
        # Slice top max_items and pad to ensure exactly max_items length
        if axis == 2:
            padded_coll = ak.pad_none(sorted_coll[:, :, :max_items], max_items, axis=axis)
        else:
            padded_coll = ak.pad_none(sorted_coll[:, :max_items], max_items, axis=axis)
        
        # Flatten into new columns
        for i in range(max_items):
            if axis == 2:
                item_i = padded_coll[:, :, i]
            else:
                item_i = padded_coll[:, i]
            for field in item_i.fields:
                col_name = f"{coll}_{i}_{field}"
                # Replace None with 0 for each field
                new_cols[col_name] = ak.fill_none(item_i[field], 0)
                
    # Zip all new columns simultaneously instead of assigning one by one (huge speedup)
    return ak.zip(new_cols, depth_limit=pt.ndim-1)

def flatten_tkel_collection(data):
    """
    Expands the nested 'tkel' collection into individual top-level fields
    and removes the original 'tkel' collection.
    """
    new_cols = {}
    for field in data.fields:
        if field == "tkel":
            for tkel_field in data.tkel.fields:
                new_cols[f"tkel_{tkel_field}"] = data.tkel[tkel_field]
        else:
            new_cols[field] = data[field]
            
    # determining depth_limit
    ndim = data.tkel.pt.ndim
    return ak.zip(new_cols, depth_limit=ndim)

def pair_tkel_puppi_and_filter(tkel, puppi, dr_cut=0.5):
    """
    Given tkel and puppi Awkward arrays, creates a cartesian product pairing 
    every TkEl with every L1PuppiCand. Calculates deta, dphi, and dR,
    adds them to the puppi records, and filters puppi candidates within dr_cut.
    """
    pairs = ak.cartesian({"tkel": tkel, "puppi": puppi}, nested=True)
    deta = pairs.tkel.eta - pairs.puppi.eta
    dphi = pairs.tkel.phi - pairs.puppi.phi
    dphi = dphi - 2 * np.pi * np.round(dphi / (2 * np.pi))
    dr = np.sqrt(deta**2 + dphi**2)
    
    # Add dR to the puppi records in pairs
    pairs["puppi", "dR"] = dr
    pairs["puppi", "deta"] = deta
    pairs["puppi", "dphi"] = dphi
    
    # Filter based on dR
    mask = (dr > 0) & (dr < dr_cut)
    filtered_puppi = pairs.puppi[mask]
    
    # Return tkel paired with matched mpuppi
    return ak.zip(
        {"tkel": tkel, "mpuppi": filtered_puppi},
        depth_limit = tkel.pt.ndim
    )
