import awkward as ak
import numpy as np
import uproot


def flatten_per_tkel(filename: str, *,
                     tkeltag = "Pt5_EB_",
                     puppitag = "Pt1_TkEleL2Pt5EB_0p0dR0p5_"):

    try:
        file = uproot.open(filename)['snapshot']
    except:
        print(f'Unable to snapshot tree in file: {filename}')
    print(f'Opened {filename}')

    # List the exact branches we need to save memory and time
    tkel_branches = {
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
    
    puppi_branches = {
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
    
    branches_to_load = list(tkel_branches.values()) + list(puppi_branches.values())
    data = file.arrays(filter_name=branches_to_load)
    print('Converted data to awkward arrays')
    
    tkel = ak.zip({k: data[v] for k, v in tkel_branches.items()})
    puppi = ak.zip({k: data[v] for k, v in puppi_branches.items()})

    # Create a cartesian product to pair every TkEl with every L1PuppCand
    pairs = ak.cartesian({"tkel": tkel, "puppi": puppi}, nested=True)
    deta = pairs.tkel.eta - pairs.puppi.eta
    dphi = pairs.tkel.phi - pairs.puppi.phi
    dphi = dphi - 2 * np.pi * np.round(dphi / (2 * np.pi))
    dr = np.sqrt(deta**2 + dphi**2)

    # Add dR to the puppi records in pairs
    pairs["puppi", "dR"] = dr
    pairs["puppi", "deta"] = deta
    pairs["puppi", "dphi"] = dphi

    # Filter based on dr
    mask = (dr > 0) & (dr < 0.5)
    filtered_puppi = pairs.puppi[mask]

    # Reassociate TkEl with matched L1PuppiCand
    tkel_with_puppi = ak.zip(
        {"tkel": tkel, "mpuppi": filtered_puppi},
        depth_limit = 2
    )

    per_tkel_data = ak.flatten(tkel_with_puppi, axis=1)
    return per_tkel_data


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
            
        # Sort by pT descending
        sorted_coll = data[coll][ak.argsort(data[coll].pt, ascending=False)]
        
        # Slice top max_items and pad to ensure exactly max_items length
        padded_coll = ak.pad_none(sorted_coll[:, :max_items], max_items, axis=1)
        
        # Flatten into new columns
        for i in range(max_items):
            item_i = padded_coll[:, i]
            for field in item_i.fields:
                col_name = f"{coll}_{i}_{field}"
                # Replace None with 0 for each field
                new_cols[col_name] = ak.fill_none(item_i[field], 0)
                
    # Zip all new columns simultaneously instead of assigning one by one (huge speedup)
    return ak.zip(new_cols, depth_limit=1)


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
            
    return ak.zip(new_cols, depth_limit=1)


def prepare_ml_data(signal_file: str, bkg_file: str):
    print("Processing signal...")
    sig_data = flatten_per_tkel(signal_file,
                                tkeltag="Pt5_EB_MCH_",
                                puppitag="Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_")
    sig_data = split_puppi_by_pdgid(sig_data)
    sig_data = flatten_puppi_collections(sig_data, max_items=5)
    sig_data = flatten_tkel_collection(sig_data)
    sig_data["label"] = 1
    print(len(sig_data))
    
    print("Processing background...")
    bkg_data = flatten_per_tkel(bkg_file)
    bkg_data = split_puppi_by_pdgid(bkg_data)
    bkg_data = flatten_puppi_collections(bkg_data, max_items=5)
    bkg_data = flatten_tkel_collection(bkg_data)
    bkg_data["label"] = 0
    print(len(bkg_data))
    
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
    train_data, test_data = prepare_ml_data(
        '../../DY_PU200_EB_snapshot.root',
        '../../MinBias_EB_snapshot.root'
    )

    print(len(train_data))
    print(len(test_data))
