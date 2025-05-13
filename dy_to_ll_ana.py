from ROOT import RDataFrame

from calc_puppi_iso import recalculate_puppi_iso
import rdf_generic as rdf_g
import utilities as ut
from utilities import sufEl, sufGen, sufPu


@ut.time_eval
def dy_to_ll_ana_main(df:RDataFrame):

    histograms = []
    
    df = df.Define(sufGen+'_n', sufGen+'_pt.size()')
    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    rdf_g.add_hists_multiplecolls(df, histograms, [sufGen, sufEl])

    df = rdf_g.define_newcollection(df, sufGen, f'abs({sufGen}_eta)<1.47 && {sufGen}_prompt==2', 'DYEB')
    rdf_g.add_hists_singlecollection(df, histograms, f'{sufGen}_DYEB')

    dfgenEB = df.Filter(f'{sufGen}_DYEB_n > 0', 'genDYEB')
    rdf_g.add_hists_multiplecolls(dfgenEB, histograms, [f'{sufGen}_DYEB', sufEl])

    #### GEN MATCH BLOCK ###
    
    dfgenEB = ut.angdiff_hists(dfgenEB, f'{sufGen}_DYEB', sufEl)
    rdf_g.add_hists_singlecollection(dfgenEB, histograms, f'{sufGen}_DYEB_{sufEl}')

    dfgenEB = ut.do_gen_match(dfgenEB,  f'{sufGen}_DYEB', sufEl)
    dfgenEB = rdf_g.define_newcollection(dfgenEB, f'{sufGen}_DYEB', f'{sufGen}_DYEB_recoidx != -1', 'MCH')
    dfgenEB = rdf_g.define_newcollection(dfgenEB, sufEl, f'{sufEl}_genidx != -1', 'MCH')
    rdf_g.add_hists_multiplecolls(dfgenEB, histograms, [f'{sufGen}_DYEB_MCH', f'{sufEl}_MCH'])

    dfgenEB = ut.angdiff_hists(dfgenEB, f'{sufGen}_DYEB_MCH', f'{sufEl}_MCH')
    rdf_g.add_hists_multiplecolls(dfgenEB, histograms, [f'{sufGen}_DYEB_MCH_{sufEl}_MCH', sufPu])

    ### END OF GEN MATCH BLOCK ###

    # Calculate puppi iso
    dfgenEB = recalculate_puppi_iso(dfgenEB, f'{sufEl}_MCH', sufPu)
    rdf_g.add_hists_singlecollection(dfgenEB, histograms, f'{sufEl}_MCH_Re', 'dRmin\\d_\\d{1,2}')

    # dfgenEB.Describe().Print()

    return histograms
