from an_specific_utilities import SampleRDFManager
import my_py_generic_utils as ut
import rdf_generic as rdf_g

@ut.time_eval
def post_analysis_persample(anamanager: SampleRDFManager):
    sname = anamanager.s_name

    if sname == 'DY_PU200':
        print(f"\nRunning post-analysis for {sname} sample...")

        DYPEB_df = anamanager.get_dataframe('DYPEB')
        # cols = ['DYPEB_TkEleL2_EB_MCH_pt', 'DYPEB_TkEleL2_EB_MCH_tkIso', 'DYPEE_TkEleL2_EE_MCH_pt', 'DYPEE_TkEleL2_EE_MCH_tkIso']
        cols = ['TkEleL2_EB_MCH_pt', 'TkEleL2_EB_MCH_tkIso']
        rdf_g.save_rdf_snapshot_to_parquet(DYPEB_df, cols, sname)