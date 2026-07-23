import awkward as ak
import numpy as np
import uproot

def flatten_per_tkel(filename: str):

    try:
        file = uproot.open(filename)['snapshot']
    except:
        print(f'Unable to snapshot tree in file: {filename}')
    print(f'Opened {filename}')

    data = file.arrays()
    print('Converted data to awkward arrays')
    tkel = ak.zip({
        "pt": data["TkEleL2_Pt5_EB_MCH_pt"],
        "tkIso": data["TkEleL2_Pt5_EB_MCH_tkIso"],
        "absTkIso": data["TkEleL2_Pt5_EB_MCH_absTkIso"],
        "puppiIso": data["TkEleL2_Pt5_EB_MCH_puppiIso"],
        "reisotot_dRmin0_03_puppiIso": data["TkEleL2_Pt5_EB_MCH_reisotot_dRmin0_03_puppiIso"],
        "reisochg_dRmin0_03_puppiIso": data["TkEleL2_Pt5_EB_MCH_reisochg_dRmin0_03_puppiIso"],
        "eta": data["TkEleL2_Pt5_EB_MCH_eta"],
        "phi": data["TkEleL2_Pt5_EB_MCH_phi"],
        "caloEta": data["TkEleL2_Pt5_EB_MCH_caloEta"],
        "caloPhi": data["TkEleL2_Pt5_EB_MCH_caloPhi"],
        "tkPt": data["TkEleL2_Pt5_EB_MCH_tkPt"],
        "tkEta": data["TkEleL2_Pt5_EB_MCH_tkEta"],
        "tkPhi": data["TkEleL2_Pt5_EB_MCH_tkPhi"],
        "charge": data["TkEleL2_Pt5_EB_MCH_charge"],
        "hwQual": data["TkEleL2_Pt5_EB_MCH_hwQual"],
        "vz": data["TkEleL2_Pt5_EB_MCH_vz"]
    })
    puppi = ak.zip({
        "pt": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_pt"],
        "eta": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_eta"],
        "phi": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_phi"],
        "mass": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_mass"],
        "charge": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_charge"],
        "dxy": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_dxy"],
        "hwDxy": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_hwDxy"],
        "hwTkQuality": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_hwTkQuality"],
        "pdgId": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_pdgId"],
        "puppiWeight": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_puppiWeight"],
        "z0": data["L1PuppiCands_Pt1_TkEleL2Pt5EBMCH_0p0dR0p5_z0"]
    })

    # Create a cartesian product to pair every TkEl with every L1PuppCand
    pairs = ak.cartesian({"tkel": tkel, "puppi": puppi}, nested=True)
    deta = pairs.tkel.eta - pairs.puppi.eta
    dphi = pairs.tkel.phi - pairs.puppi.phi
    dphi = dphi - 2 * np.pi * np.round(dphi / (2 * np.pi))
    dr = np.sqrt(deta**2 + dphi**2)

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

d = flatten_per_tkel('../../DY_PU200_EB_snapshot.root')