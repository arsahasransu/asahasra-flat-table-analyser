from an_specific_utilities import SampleRDFManager
import pypkg.my_py_generic_utils as ut
import rdf_generic as rdf_g

@ut.time_eval
def post_analysis_persample(anamanager: SampleRDFManager):
    sname = anamanager.s_name

    if sname == 'DY_PU200':
        print(f"\nRunning post-analysis for {sname} sample...")

        # STEP 3_0_0, X_X_X: For ROC curves
        DYPEB_df = anamanager.get_dataframe('DYPMEBP')
        cols = ['TkEleL2_Pt5_EB_MCH_pt', 'TkEleL2_Pt5_EB_MCH_tkIso', 'TkEleL2_Pt5_EB_MCH_absTkIso',
                'TkEleL2_Pt5_EB_MCH_puppiIso',
                'TkEleL2_Pt5_EB_MCH_reisotot_dRmin0_03_puppiIso',
                'TkEleL2_Pt5_EB_MCH_reisochg_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EB_MCH_reisonut_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EB_MCH_reisotot_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EB_MCH_reisochg_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EB_MCH_reisonut_dRmin0_01_puppiIso',
            ]
        cols += add_ntuples_for_tkel('EB_MCH')
        cols += add_ntuples_for_puppi('EBMCH')
        rdf_g.save_rdf_snapshot(DYPEB_df, cols, f"{sname}_EB", recreate=True)

        # STEP 3_0_0, X_X_X: For ROC curves
        DYPEE_df = anamanager.get_dataframe('DYPMEEP')
        cols = ['TkEleL2_Pt5_EE_MCH_pt', 'TkEleL2_Pt5_EE_MCH_tkIso', 'TkEleL2_Pt5_EE_MCH_absTkIso',
                'TkEleL2_Pt5_EE_MCH_puppiIso',
                'TkEleL2_Pt5_EE_MCH_reisotot_dRmin0_03_puppiIso',
                'TkEleL2_Pt5_EE_MCH_reisochg_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EE_MCH_reisonut_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EE_MCH_reisotot_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EE_MCH_reisochg_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EE_MCH_reisonut_dRmin0_01_puppiIso',
            ]  
        cols += add_ntuples_for_tkel('EE_MCH')
        cols += add_ntuples_for_puppi('EEMCH')
        rdf_g.save_rdf_snapshot(DYPEE_df, cols, f"{sname}_EE", recreate=True)

    if sname == "MinBias":
        print(f"\nRunning post-analysis for {sname} sample...")

        # STEP 3_0_0, X_X_X: For ROC curves
        EB_df = anamanager.get_dataframe('EBP')
        cols = ['TkEleL2_Pt5_EB_pt', 'TkEleL2_Pt5_EB_tkIso', 'TkEleL2_Pt5_EB_absTkIso',
                'TkEleL2_Pt5_EB_puppiIso',
                'TkEleL2_Pt5_EB_reisotot_dRmin0_03_puppiIso',
                'TkEleL2_Pt5_EB_reisochg_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EB_reisonut_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EB_reisotot_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EB_reisochg_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EB_reisonut_dRmin0_01_puppiIso',
            ]
        cols += add_ntuples_for_tkel('EB')
        cols += add_ntuples_for_puppi('EB')
        rdf_g.save_rdf_snapshot(EB_df, cols, f"{sname}_EB", recreate=True)

        # STEP 3_0_0, X_X_X: For ROC curves
        EE_df = anamanager.get_dataframe('EEP')
        cols = ['TkEleL2_Pt5_EE_pt', 'TkEleL2_Pt5_EE_tkIso', 'TkEleL2_Pt5_EE_absTkIso',
                'TkEleL2_Pt5_EE_puppiIso',
                'TkEleL2_Pt5_EE_reisotot_dRmin0_03_puppiIso',
                'TkEleL2_Pt5_EE_reisochg_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EE_reisonut_dRmin0_03_puppiIso',
                # 'TkEleL2_Pt5_EE_reisotot_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EE_reisochg_dRmin0_01_puppiIso',
                # 'TkEleL2_Pt5_EE_reisonut_dRmin0_01_puppiIso',
        ]
        cols += add_ntuples_for_tkel('EE')
        cols += add_ntuples_for_puppi('EE')
        rdf_g.save_rdf_snapshot(EE_df, cols, f"{sname}_EE", recreate=True)


def add_ntuples_for_tkel(eregion):

    cols = []
    vars = ['eta', 'phi', 'caloEta', 'caloPhi', 'tkPt', 'tkEta', 'tkPhi', 'charge', 'hwQual', 'vz']

    for var in vars:
        collname = f'TkEleL2_Pt5_{eregion}_{var}'
        cols.append(collname)

    return cols


def add_ntuples_for_puppi(eregion):
    cols = []
    collnamepart = f'L1PuppiCands_Pt1_TkEleL2Pt5{eregion}_0p0dR0p5'

    # vars = ['pt', 'eta', 'phi']
    vars = ['pt', 'eta', 'phi', 'mass', 'charge', 'dxy', 'hwDxy', 'hwTkQuality', 'pdgId', 'puppiWeight', 'z0']

    for var in vars:
        collname = f'{collnamepart}_{var}'
        cols.append(collname)

    print(cols)
    return cols
