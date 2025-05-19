from ROOT import RDataFrame

import calc_puppi_iso as reiso
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

    # Checked for pT sorting of the TkEl collection
    # dfgenEB = dfgenEB.Define(f'{sufEl}_MCH_isptsorted', f'checksorting<float>({sufEl}_MCH_pt, true)')
    # histograms.append( dfgenEB.Histo1D((f'genDYEB_{sufEl}_MCH_isptsorted', 'ptsorting', 6, -2, 4),\
    # f'{sufEl}_MCH_isptsorted') )

    # Calculate puppi iso
    # dfgenEB = reiso.recalculate_puppi_iso(dfgenEB, f'{sufEl}_MCH', sufPu)
    # rdf_g.add_hists_singlecollection(dfgenEB, histograms, f'{sufEl}_MCH_Re', 'dRmin\\d_\\d{1,2}')

    # Separate the lead and sub-lead gen matched TkEl
    dfgEBel1 = dfgenEB.Filter(f'{sufEl}_MCH_n >= 1', 'genDYEBTkEl1')
    dfgEBel1 = rdf_g.define_newcollection(dfgEBel1, f'{sufEl}_MCH', '0', 'El1')
    rdf_g.add_hists_singlecollection(dfgEBel1, histograms, f'{sufEl}_MCH_El1')
    # Puppi isolation and components by pdgId
    dfgEBel1 = reiso.recalcpuppiiso_comps_oneel(dfgEBel1, f'{sufEl}_MCH_El1', sufPu)
    rdf_g.add_hists_singlecollection(dfgEBel1, histograms, f'{sufEl}_MCH_El1_Re',\
                                     'dRmin\\d_\\d{1,2}_dRmax\\d_\\d{1,2}_[a-z0-9]+')

    # dfgEBel2 = dfgenEB.Filter(f'{sufEl}_MCH_n >= 2', 'genDYEBTkEl2')
    # dfgEBel2 = rdf_g.define_newcollection(dfgEBel2, f'{sufEl}_MCH', '1', 'El2')
    # rdf_g.add_hists_singlecollection(dfgEBel2, histograms, f'{sufEl}_MCH_El2')

    # Calculate the invariant mass
    # dfgEBel2 = dfgEBel2.Define(f'{sufEl}_MCH_m', f"ROOT::VecOps::RVec<float>({sufEl}_MCH_pt.size(), 0.000511)")
    # dfgEBel2 = dfgEBel2.Define(f'{sufEl}_MCH_Mee', f'ROOT::VecOps::InvariantMass<float>({sufEl}_MCH_pt,\
    # {sufEl}_MCH_eta, {sufEl}_MCH_phi, {sufEl}_MCH_m)')
    # histograms.append( dfgEBel2.Histo1D((f'genDYEBTkEl2_{sufEl}_MCH_Mee', 'mass [GeV]', 201, -1, 200),\
    # f'{sufEl}_MCH_Mee') )

    # dfgenEB_tkel2.Describe().Print()

    return histograms
