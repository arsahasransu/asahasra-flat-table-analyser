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

    data = file.arrays()
    print('Converted data to awkward arrays')
    tkel = ak.zip({
        "pt": data[f"TkEleL2_{tkeltag}pt"],
        "eta": data[f"TkEleL2_{tkeltag}eta"],
        "phi": data[f"TkEleL2_{tkeltag}phi"],
        "caloEta": data[f"TkEleL2_{tkeltag}caloEta"],
        "caloPhi": data[f"TkEleL2_{tkeltag}caloPhi"],
        "tkPt": data[f"TkEleL2_{tkeltag}tkPt"],
        "tkEta": data[f"TkEleL2_{tkeltag}tkEta"],
        "tkPhi": data[f"TkEleL2_{tkeltag}tkPhi"],
        "charge": data[f"TkEleL2_{tkeltag}charge"],
        "hwQual": data[f"TkEleL2_{tkeltag}hwQual"],
        "vz": data[f"TkEleL2_{tkeltag}vz"]
    })
    puppi = ak.zip({
        "pt": data[f"L1PuppiCands_{puppitag}pt"],
        "eta": data[f"L1PuppiCands_{puppitag}eta"],
        "phi": data[f"L1PuppiCands_{puppitag}phi"],
        "mass": data[f"L1PuppiCands_{puppitag}mass"],
        "charge": data[f"L1PuppiCands_{puppitag}charge"],
        "dxy": data[f"L1PuppiCands_{puppitag}dxy"],
        "hwDxy": data[f"L1PuppiCands_{puppitag}hwDxy"],
        "hwTkQuality": data[f"L1PuppiCands_{puppitag}hwTkQuality"],
        "pdgId": data[f"L1PuppiCands_{puppitag}pdgId"],
        "puppiWeight": data[f"L1PuppiCands_{puppitag}puppiWeight"],
        "z0": data[f"L1PuppiCands_{puppitag}z0"]
    })

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
    
    # Remove original 'mpuppi' if present to save memory
    if "mpuppi" in data.fields:
        data = data[[field for field in data.fields if field != "mpuppi"]]
        
    for coll in collections:
        if coll not in data.fields:
            print(f"Warning, did not find {field} in data")
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
                data[col_name] = ak.fill_none(item_i[field], 0)
                
        # Remove the nested collection
        data = data[[f for f in data.fields if f != coll]]
        
    return data

def prepare_ml_data(signal_file: str, bkg_file: str):
    print("Processing signal...")
    sig_data = flatten_per_tkel(signal_file,
                                tkeltag="Pt5_EB_MCH_",
                                puppitag="Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_")
    sig_data = split_puppi_by_pdgid(sig_data)
    sig_data = flatten_puppi_collections(sig_data, max_items=5)
    sig_data["label"] = 1
    print(len(sig_data))
    
    print("Processing background...")
    bkg_data = flatten_per_tkel(bkg_file)
    bkg_data = split_puppi_by_pdgid(bkg_data)
    bkg_data = flatten_puppi_collections(bkg_data, max_items=5)
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
