from ROOT import RDataFrame

from calc_puppi_iso import recalculate_puppi_iso
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
    rdf_g.add_hists_singlecollection(df, histograms, f'{sufEl}_TkElEBPt10')

    dfTkElEBPt10 = df.Filter(f'{sufEl}_TkElEBPt10_n > 0', 'tkelEBpt10')
    rdf_g.add_hists_multiplecolls(dfTkElEBPt10, histograms, [f'{sufEl}_TkElEBPt10', sufPu])

    dfTkElEBPt10 = recalculate_puppi_iso(dfTkElEBPt10, f'{sufEl}_TkElEBPt10', sufPu)
    rdf_g.add_hists_singlecollection(dfTkElEBPt10, histograms, f'{sufEl}_TkElEBPt10_Re', 'dRmin\\d_\\d{1,2}')

    return histograms