from ROOT import RDataFrame

import rdf_generic as rdf_g

def post_analysis_persample(sname: str, df: RDataFrame):
    
    df.Describe().Print()

    if sname == 'DY_PU200':
        # cols = ['DYPEB_TkEleL2_EB_MCH_pt', 'DYPEB_TkEleL2_EB_MCH_tkIso', 'DYPEE_TkEleL2_EE_MCH_pt', 'DYPEE_TkEleL2_EE_MCH_tkIso']
        cols = ['DYPEB_TkEleL2_EB_MCH_tkIso']
        rdf_g.save_rdf_snapshot_to_parquet(df, cols, sname)