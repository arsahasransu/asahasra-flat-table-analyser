from ROOT import RDataFrame

import calc_puppi_iso as reiso
import rdf_generic as rdf_g
import an_specific_utilities as ut
from an_specific_utilities import sufEl, sufGen, sufPu
from an_specific_utilities import add_puppicands_by_pdg


@ut.time_eval
def dy_to_ll_ana_main(df: RDataFrame):

    histograms = []

    df = df.Define(sufGen+'_n', sufGen+'_pt.size()')
    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    # Step 1_0_0: Enable for plots in "Gen properties" section
    ##########################################################
    # rdf_g.add_hists_multiplecolls(df, histograms, [sufGen, sufEl, sufPu])

    # df_genstudy1 = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==1', 'GENST1')
    # rdf_g.add_hists_singlecollection(df_genstudy1, histograms, f'{sufGen}_GENST1')

    # df_genstudy0 = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==0', 'GENST0')
    # rdf_g.add_hists_singlecollection(df_genstudy0, histograms, f'{sufGen}_GENST0')
    ##########################################################

    df = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt>=2 && abs({sufGen}_eta)<2.6', 'DYP')

    dfGenP = df.Filter(f'{sufGen}_DYP_n > 0 && {sufEl}_n > 0', 'genDYP')
    # Step 1_2_0: Enable for plots in "Gen selection" section
    ##########################################################
    # rdf_g.add_hists_multiplecolls(dfGenP, histograms, [f'{sufGen}_DYP', sufEl, sufPu])
    # add_puppicands_by_pdg(dfGenP, histograms, '')
    ##########################################################

    # SPLIT ETA REGIONS BASED ON TKEL
    dfGenP = rdf_g.define_newcollection(dfGenP, sufEl, f'abs({sufEl}_eta) <= 1.4', 'EB')
    dfGenP = rdf_g.define_newcollection(dfGenP, sufEl, f'abs({sufEl}_eta) > 1.4 && abs({sufEl}_eta) <= 1.6', 'EM')
    dfGenP = rdf_g.define_newcollection(dfGenP, sufEl, f'abs({sufEl}_eta) > 1.6 && abs({sufEl}_eta) <= 2.1', 'EE')
    dfGenP = rdf_g.define_newcollection(dfGenP, sufEl, f'abs({sufEl}_eta) > 2.1', 'EF')

    dRcuts = {'EB': 0.02, 'EM': 0.025, 'EE': 0.04, 'EF': 0.05}

    for ERegion in ['EB', 'EM', 'EE', 'EF']:
        sufGenDYP = f'{sufGen}_DYP'
        sufElER = f'{sufEl}_{ERegion}'
        dfGenER = dfGenP.Filter(f'{sufElER}_n > 0', f'genDYPel{ERegion}')

        # Step 2_0_0: GEN MATCH BLOCK
        ##########################################################
        # PRE GEN MATCH
        # dfGenER = ut.angdiff_hists(dfGenER, sufGenDYP, sufElER)
        # rdf_g.add_hists_singlecollection(dfGenER, histograms, f'{sufGenDYP}_{sufElER}')
        # ut.add_genmatching_efficiency_with_dRcut(histograms, f'genDYPel{ERegion}_{sufGenDYP}_{sufElER}')

        # GEN MATCH
        dfGenER = ut.do_gen_match(dfGenER, sufGenDYP, sufElER, dRcuts[ERegion])
        dfGenER = rdf_g.define_newcollection(dfGenER, sufGenDYP, f'{sufGenDYP}_recoidx != -1', 'MCH')
        dfGenER = rdf_g.define_newcollection(dfGenER, sufElER, f'{sufElER}_genidx != -1', 'MCH')\
        
        sufGenMch = f'{sufGenDYP}_MCH'
        sufElMch = f'{sufElER}_MCH'
        rdf_g.add_hists_multiplecolls(dfGenER, histograms, [sufGenMch, sufElMch])

        # POST GEN MATCH
        # dfGenER = ut.angdiff_hists(dfGenER, sufGenMch, sufElMch)
        # rdf_g.add_hists_singlecollection(dfGenER, histograms, f'{sufGenMch}_{sufElMch}')
        ##########################################################

    return histograms

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
    rdf_g.add_hists_singlecollection(dfgEBel1, histograms, f'{sufEl}_MCH_El1_Re',
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
