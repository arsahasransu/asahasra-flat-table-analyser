from ROOT import RDataFrame

import an_specific_utilities as anut
from an_specific_utilities import sufEl, sufGen, sufPu
from an_specific_utilities import add_puppicands_by_pdg
import calc_puppi_iso as reiso
import my_py_generic_utils as ut
import rdf_generic as rdf_g


@ut.time_eval
def dy_to_ll_ana_main(df: RDataFrame):

    histograms = []

    df = df.Define(sufGen+'_n', sufGen+'_pt.size()')
    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    # STEP 1_0_0: Enable for plots in "Gen properties" section
    ##########################################################
    # rdf_g.add_hists_multiplecolls(df, histograms, [sufGen, sufEl, sufPu])

    # df_genstudy2 = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==2', 'GENST2')
    # rdf_g.add_hists_singlecollection(df_genstudy2, histograms, f'{sufGen}_GENST2')

    # df_genstudy1 = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==1', 'GENST1')
    # rdf_g.add_hists_singlecollection(df_genstudy1, histograms, f'{sufGen}_GENST1')

    # df_genstudy0 = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==0', 'GENST0')
    # rdf_g.add_hists_singlecollection(df_genstudy0, histograms, f'{sufGen}_GENST0')
    
    df = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==2 && abs({sufGen}_eta)<2.5', 'DYP')
    # rdf_g.add_hists_singlecollection(df, histograms, f'{sufGen}_DYP')
    ##########################################################

    dfGenPS = df.Filter(f'{sufGen}_DYP_n > 0', 'genDYP')
    ut.create_rdf_checkpint(df, dfGenPS, "Applying selection: >0 Gen electron with prompt status 2 and |eta| < 2.5...")

    dfGenP = dfGenPS.Filter(f'{sufGen}_DYP_n > 0 && {sufEl}_n > 0', 'genDYP')
    ut.create_rdf_checkpint(dfGenPS, dfGenP, "Applying selection: >0 reconstructed TkEl...")

    # STEP 1_2_0: Enable for plots in "Gen selection" section
    ##########################################################
    # rdf_g.add_hists_singlecollection(dfGenP, histograms, sufEl)
    # add_puppicands_by_pdg(dfGenP, histograms, '')
    ##########################################################

    # SPLIT ETA REGIONS BASED ON TKEL
    dfGenP = rdf_g.define_newcollection(dfGenP, sufEl, f'abs({sufEl}_eta) <= 1.5', 'EB')
    dfGenP = rdf_g.define_newcollection(dfGenP, sufEl, f'abs({sufEl}_eta) > 1.5 && abs({sufEl}_eta) <= 2.5', 'EE')

    gen_dRcuts = {'EB': 0.04, 'EE': 0.05}

    for ERegion in ['EB', 'EE']:
    # for ERegion in ['EB']:
        sufGenDYP = f'{sufGen}_DYP'
        sufElER = f'{sufEl}_{ERegion}'
        dfGenER = dfGenP.Filter(f'{sufElER}_n > 0', f'DYP{ERegion}')
        ut.create_rdf_checkpint(dfGenP, dfGenER, f"Applying selection: > 0 TkEl in region {ERegion}...")

        # STEP 2_0_0: GEN MATCH BLOCK
        #########################################################
        # PRE GEN MATCH
        dfGenER = anut.angdiff_hists(dfGenER, sufGenDYP, sufElER)
        rdf_g.add_hists_singlecollection(dfGenER, histograms, f'{sufGenDYP}_{sufElER}')
        anut.add_genmatching_efficiency_with_dRcut(histograms, f'DYP{ERegion}_{sufGenDYP}_{sufElER}')

        # GEN MATCH
        dfGenER = anut.do_gen_match(dfGenER, sufGenDYP, sufElER, gen_dRcuts[ERegion])
        dfGenER = rdf_g.define_newcollection(dfGenER, sufGenDYP, f'{sufGenDYP}_recoidx != -1', 'MCH')
        dfGenER = rdf_g.define_newcollection(dfGenER, sufElER, f'{sufElER}_genidx != -1', 'MCH')\
        
        sufGenMch = f'{sufGenDYP}_MCH'
        sufElMch = f'{sufElER}_MCH'

        # POST GEN MATCH
        dfGenER = anut.angdiff_hists(dfGenER, sufGenMch, sufElMch)
        rdf_g.add_hists_singlecollection(dfGenER, histograms, f'{sufGenMch}_{sufElMch}')
        ##########################################################

        # Filter for atleast one gen-match TkEl in the defined eta region
        # dfGenER = dfGenER.Filter(f'{sufElMch}_n > 0', f'genDYPTkElMCH{ERegion}')
        # rdf_g.add_hists_multiplecolls(dfGenER, histograms, [sufGenMch, sufElMch])
        # rdf_g.add_hists_multiplecolls(dfGenER, histograms, [sufElMch])

        # # # Checked for pT sorting of the TkEl collection
        # dfGenER = dfGenER.Define(f'{sufElMch}_pt_sorted', f'checksorting<float>({sufElMch}_pt, true)')
        # histograms.append( dfGenER.Histo1D((f'{sufElMch}_pt_sorted', 'pt_sort', 6, -2, 4), f'{sufElMch}_pt_sorted') )
        # cp_pt_sort = anut.check_histogram_for_value(histograms, f'{sufElMch}_pt_sorted', bin=3) # bin 3 = value 0
        # if cp_pt_sort:
        #     raise RuntimeError("The Track Electrons pT's are not sorted from highest to lowest")

        # For each TkEl, observe the feature of PuppiCandidates around the TkEl in annular dR window
        # add_puppicands_by_pdg(dfGenER, histograms, '', tkelobj=sufElMch)
        # dfGenER = anut.make_puppi_by_angdiff_from_tkel(dfGenER, sufElMch, histograms)

        # dfGenER.Describe().Print()

    return histograms

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
