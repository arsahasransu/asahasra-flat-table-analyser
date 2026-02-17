from an_specific_utilities import SampleRDFManager
import my_py_generic_utils as ut
import rdf_generic as rdf_g

@ut.time_eval
def post_analysis_persample(anamanager: SampleRDFManager):
    sname = anamanager.s_name

    if sname == 'DY_PU200':
        print(f"\nRunning post-analysis for {sname} sample...")

        # STEP 3_0_0, X_X_X: For ROC curves
        DYPEB_df = anamanager.get_dataframe('DYPMEB')
        cols = ['TkEleL2_EB_MCH_pt', 'TkEleL2_EB_MCH_tkIso', 'TkEleL2_EB_MCH_puppiIso',
                'TkEleL2_EB_MCH_reisotot_dRmin0_01_puppiIso',
                'TkEleL2_EB_MCH_reisotot_dRmin0_01_absPuppiIso']
        rdf_g.save_rdf_snapshot_to_pkl(DYPEB_df, cols, sname, recreate=True)

        # STEP 3_0_0, X_X_X: For ROC curves
        DYPEE_df = anamanager.get_dataframe('DYPMEE')
        cols = ['TkEleL2_EE_MCH_pt', 'TkEleL2_EE_MCH_tkIso', 'TkEleL2_EE_MCH_puppiIso',
                'TkEleL2_EE_MCH_reisotot_dRmin0_01_puppiIso',
                'TkEleL2_EE_MCH_reisotot_dRmin0_01_absPuppiIso']
        rdf_g.save_rdf_snapshot_to_pkl(DYPEE_df, cols, sname)

    if sname == "MinBias":
        print(f"\nRunning post-analysis for {sname} sample...")

        # STEP 3_0_0, X_X_X: For ROC curves
        EB_df = anamanager.get_dataframe('EB')
        cols = ['TkEleL2_EB_pt', 'TkEleL2_EB_tkIso', 'TkEleL2_EB_puppiIso',
                'TkEleL2_EB_reisotot_dRmin0_01_puppiIso',
                'TkEleL2_EB_reisotot_dRmin0_01_absPuppiIso']
        rdf_g.save_rdf_snapshot_to_pkl(EB_df, cols, sname, recreate=True)

        # STEP 3_0_0, X_X_X: For ROC curves
        EE_df = anamanager.get_dataframe('EE')
        cols = ['TkEleL2_EE_pt', 'TkEleL2_EE_tkIso', 'TkEleL2_EE_puppiIso',
                'TkEleL2_EE_reisotot_dRmin0_01_puppiIso',
                'TkEleL2_EE_reisotot_dRmin0_01_absPuppiIso']
        rdf_g.save_rdf_snapshot_to_pkl(EE_df, cols, sname)