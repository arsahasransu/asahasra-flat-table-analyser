from ROOT import RDataFrame

import rdf_generic as rdf_g
import utilities as ut
from utilities import sufEl, sufPu


@ut.time_eval
def qcd_ana_main(df:RDataFrame):

    histograms = []
    
    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    rdf_g.add_hists_multiplecolls(df, histograms, [sufEl, sufPu])

    df = rdf_g.define_newcollection(df, sufEl, f'{sufEl}_pt > 10 && abs({sufEl}_eta) < 1.479', 'TkElEBPt10')
    rdf_g.add_hists_singlecollection(df, histograms, sufEl, 'TkElEBPt10')

    dfTkElEBPt10 = df.Filter(f'{sufEl}_TkElEBPt10_n > 0', 'tkelEBpt10')
    rdf_g.add_hists_multiplecolls(dfTkElEBPt10, histograms, [f'{sufEl}_TkElEBPt10', sufPu])

    # df_TkElEB = get_iso(df_TkElEB, f'{sufEl}TkElPt10EB', sufPu)
    # dRmin_list = np.arange(0.01, 0.2, 0.01)
    # for dRmin in dRmin_list:
    #     histograms.append(df_TkElEB.Histo1D((f'{sufEl}TkElPt10EB_recalcPuppiIso_dRmin{str(dRmin).replace('.', '_')}', 
    #                                          'recalcPuppiIso', 10000, -2, 8), f'{sufEl}TkElPt10EB_recalcPuppiIso_dRmin{str(dRmin).replace('.', '_')}'))

    return histograms