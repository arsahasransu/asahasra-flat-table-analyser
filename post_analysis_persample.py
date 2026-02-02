from an_specific_utilities import SampleRDFManager
import my_py_generic_utils as ut
import rdf_generic as rdf_g

@ut.time_eval
def post_analysis_persample(anamanager: SampleRDFManager):
    sname = anamanager.s_name

    if sname == 'DY_PU200':
        print(f"\nRunning post-analysis for {sname} sample...")

        DYPEB_df = anamanager.get_dataframe('DYPEB')
        cols = ['TkEleL2_EB_MCH_pt', 'TkEleL2_EB_MCH_tkIso', 'TkEleL2_EB_MCH_puppiIso']
        rdf_g.save_rdf_snapshot_to_pkl(DYPEB_df, cols, sname, recreate=True)

        DYPEE_df = anamanager.get_dataframe('DYPEE')
        cols = ['TkEleL2_EE_MCH_pt', 'TkEleL2_EE_MCH_tkIso', 'TkEleL2_EE_MCH_puppiIso']
        rdf_g.save_rdf_snapshot_to_pkl(DYPEE_df, cols, sname)

    if sname == "MinBias":
        print(f"\nRunning post-analysis for {sname} sample...")

        EB_df = anamanager.get_dataframe('EB')
        cols = ['TkEleL2_EB_pt', 'TkEleL2_EB_tkIso', 'TkEleL2_EB_puppiIso']
        rdf_g.save_rdf_snapshot_to_pkl(EB_df, cols, sname, recreate=True)

        EE_df = anamanager.get_dataframe('EE')
        cols = ['TkEleL2_EE_pt', 'TkEleL2_EE_tkIso', 'TkEleL2_EE_puppiIso']
        rdf_g.save_rdf_snapshot_to_pkl(EE_df, cols, sname)