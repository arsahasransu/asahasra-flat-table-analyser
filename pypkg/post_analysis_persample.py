from an_specific_utilities import SampleRDFManager
import pypkg.my_py_generic_utils as ut
import rdf_generic as rdf_g

@ut.time_eval
def post_analysis_persample(anamanager: SampleRDFManager):
    sname = anamanager.s_name

    if sname == 'DY_PU200':
        print(f"\nRunning post-analysis for {sname} sample...")

        # STEP 3_0_0, X_X_X: For ROC curves
        DYPEB_df = anamanager.get_dataframe('DYPMEB')
        cols = ['TkEleL2_Pt5_EB_MCH_pt', 'TkEleL2_Pt5_EB_MCH_tkIso', 'TkEleL2_Pt5_EB_MCH_absTkIso',
                'TkEleL2_Pt5_EB_MCH_puppiIso',
                'TkEleL2_Pt5_EB_MCH_reisotot_dRmin0_03_puppiIso',
                'TkEleL2_Pt5_EB_MCH_reisochg_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EB_MCH_reisonut_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EB_MCH_reisotot_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EB_MCH_reisochg_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EB_MCH_reisonut_dRmin0_01_puppiIso',
            ]
        cols += add_ntuples_for_aialgos_dy_eb()
        cols += add_ntuples_for_aialgos_puppi()
        rdf_g.save_rdf_snapshot(DYPEB_df, cols, f"{sname}_EB", recreate=True)

        # STEP 3_0_0, X_X_X: For ROC curves
        DYPEE_df = anamanager.get_dataframe('DYPMEE')
        cols = ['TkEleL2_Pt5_EE_MCH_pt', 'TkEleL2_Pt5_EE_MCH_tkIso', 'TkEleL2_Pt5_EE_MCH_absTkIso',
                'TkEleL2_Pt5_EE_MCH_puppiIso',
                'TkEleL2_Pt5_EE_MCH_reisotot_dRmin0_03_puppiIso',
                'TkEleL2_Pt5_EE_MCH_reisochg_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EE_MCH_reisonut_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EE_MCH_reisotot_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EE_MCH_reisochg_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EE_MCH_reisonut_dRmin0_01_puppiIso',
            ]  
        cols += add_ntuples_for_aialgos_dy_ee()
        cols += add_ntuples_for_aialgos_puppi()
        rdf_g.save_rdf_snapshot(DYPEE_df, cols, f"{sname}_EE", recreate=True)

    if sname == "MinBias":
        print(f"\nRunning post-analysis for {sname} sample...")

        # STEP 3_0_0, X_X_X: For ROC curves
        EB_df = anamanager.get_dataframe('EB')
        cols = ['TkEleL2_Pt5_EB_pt', 'TkEleL2_Pt5_EB_tkIso', 'TkEleL2_Pt5_EB_absTkIso',
                'TkEleL2_Pt5_EB_puppiIso',
                'TkEleL2_Pt5_EB_reisotot_dRmin0_03_puppiIso',
                'TkEleL2_Pt5_EB_reisochg_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EB_reisonut_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EB_reisotot_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EB_reisochg_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EB_reisonut_dRmin0_01_puppiIso',
            ]
        cols += add_ntuples_for_aialgos_minbias_eb()
        cols += add_ntuples_for_aialgos_puppi()
        rdf_g.save_rdf_snapshot(EB_df, cols, f"{sname}_EB", recreate=True)

        # STEP 3_0_0, X_X_X: For ROC curves
        EE_df = anamanager.get_dataframe('EE')
        cols = ['TkEleL2_Pt5_EE_pt', 'TkEleL2_Pt5_EE_tkIso', 'TkEleL2_Pt5_EE_absTkIso',
                'TkEleL2_Pt5_EE_puppiIso',
                'TkEleL2_Pt5_EE_reisotot_dRmin0_03_puppiIso',
                'TkEleL2_Pt5_EE_reisochg_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EE_reisonut_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EE_reisotot_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EE_reisochg_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EE_reisonut_dRmin0_01_puppiIso',
        ]
        cols += add_ntuples_for_aialgos_minbias_ee()
        cols += add_ntuples_for_aialgos_puppi()

        rdf_g.save_rdf_snapshot(EE_df, cols, f"{sname}_EE", recreate=True)


def add_ntuples_for_aialgos_dy_eb():

    cols = ['TkEleL2_Pt5_EB_MCH_eta', 'TkEleL2_Pt5_EB_MCH_phi',
            'TkEleL2_Pt5_EB_MCH_caloEta', 'TkEleL2_Pt5_EB_MCH_caloPhi',
            'TkEleL2_Pt5_EB_MCH_tkPt', 'TkEleL2_Pt5_EB_MCH_tkEta', 'TkEleL2_Pt5_EB_MCH_tkPhi',
            'TkEleL2_Pt5_EB_MCH_charge', 'TkEleL2_Pt5_EB_MCH_hwQual', 'TkEleL2_Pt5_EB_MCH_vz']

    return cols


def add_ntuples_for_aialgos_dy_ee():

    cols = ['TkEleL2_Pt5_EE_MCH_eta', 'TkEleL2_Pt5_EE_MCH_phi',
            'TkEleL2_Pt5_EE_MCH_caloEta', 'TkEleL2_Pt5_EE_MCH_caloPhi',
            'TkEleL2_Pt5_EE_MCH_tkPt', 'TkEleL2_Pt5_EE_MCH_tkEta', 'TkEleL2_Pt5_EE_MCH_tkPhi',
            'TkEleL2_Pt5_EE_MCH_charge', 'TkEleL2_Pt5_EE_MCH_hwQual', 'TkEleL2_Pt5_EE_MCH_vz']

    return cols


def add_ntuples_for_aialgos_minbias_eb():

    cols = ['TkEleL2_Pt5_EB_eta', 'TkEleL2_Pt5_EB_phi',
            'TkEleL2_Pt5_EB_caloEta', 'TkEleL2_Pt5_EB_caloPhi',
            'TkEleL2_Pt5_EB_tkPt', 'TkEleL2_Pt5_EB_tkEta', 'TkEleL2_Pt5_EB_tkPhi',
            'TkEleL2_Pt5_EB_charge', 'TkEleL2_Pt5_EB_hwQual', 'TkEleL2_Pt5_EB_vz']
    
    return cols


def add_ntuples_for_aialgos_minbias_ee():

    cols = ['TkEleL2_Pt5_EE_eta', 'TkEleL2_Pt5_EE_phi',
            'TkEleL2_Pt5_EE_caloEta', 'TkEleL2_Pt5_EE_caloPhi',
            'TkEleL2_Pt5_EE_tkPt', 'TkEleL2_Pt5_EE_tkEta', 'TkEleL2_Pt5_EE_tkPhi',
            'TkEleL2_Pt5_EE_charge', 'TkEleL2_Pt5_EE_hwQual', 'TkEleL2_Pt5_EE_vz']

    return cols

def add_ntuples_for_aialgos_puppi():

    cols = ['L1PuppiCands_Pt1_pt', 'L1PuppiCands_Pt1_eta', 'L1PuppiCands_Pt1_phi',
            'L1PuppiCands_Pt1_mass', 'L1PuppiCands_Pt1_charge', 'L1PuppiCands_Pt1_dxy',
            'L1PuppiCands_Pt1_hwDxy', 'L1PuppiCands_Pt1_hwTkQuality', 'L1PuppiCands_Pt1_pdgId',
            'L1PuppiCands_Pt1_puppiWeight', 'L1PuppiCands_Pt1_z0']

    return cols
