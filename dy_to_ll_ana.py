from ROOT import RDataFrame

import an_specific_utilities as anut
from an_specific_utilities import sufEl, sufGen, sufPu
from an_specific_utilities import add_puppicands_by_pdg
import pypkg.calc_puppi_iso as reiso
import pypkg.my_py_generic_utils as ut
import rdf_generic as rdf_g


@ut.time_eval
def dy_to_ll_ana_main(ana_man: anut.SampleRDFManager) -> anut.SampleRDFManager:

    histograms = []

    df = ana_man.parent_df
    df = df.Define(sufGen+'_n', sufGen+'_pt.size()')
    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    df = df.Define(sufEl+'_absTkIso', f'{sufEl}_pt*{sufEl}_tkIso')

    # # STEP 1_0_0: Enable for plots in "Gen properties" section
    # ##########################################################
    # rdf_g.add_hists_multiplecolls(df, histograms, [sufGen, sufEl, sufPu])

    # df_genstudy2 = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==2', 'GENST2')
    # rdf_g.add_hists_singlecollection(df_genstudy2, histograms, f'{sufGen}_GENST2')

    # df_genstudy1 = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==1', 'GENST1')
    # rdf_g.add_hists_singlecollection(df_genstudy1, histograms, f'{sufGen}_GENST1')

    # df_genstudy0 = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==0', 'GENST0')
    # rdf_g.add_hists_singlecollection(df_genstudy0, histograms, f'{sufGen}_GENST0')    
    # ##########################################################

    df = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==2 && abs({sufGen}_eta)<=2.5', 'DYP')
    df = rdf_g.define_newcollection(df, sufEl, f'{sufEl}_pt >= 15.0', 'Pt5')
    sufElPt5 = sufEl+'_Pt5'
    df = rdf_g.define_newcollection(df, sufPu, f'{sufPu}_pt >= 1.0', 'Pt1')
    sufPuPt1 = sufPu+'_Pt1'
    rdf_g.add_hists_singlecollection(df, histograms, f'{sufGen}_DYP')

    dfGenP = df.Filter(f'{sufGen}_DYP_n > 0 && {sufElPt5}_n > 0', 'genDYP')
    ana_man.add_dataframe(key='genDYP', df=dfGenP)
    ut.create_rdf_checkpint(df, dfGenP, "Applying selection: >0 Gen electron with prompt "
                                        "status 2 and |eta| < 2.5... \n>0 reconstructed TkEl...")


    # # STEP 1_2_0: Enable for plots in "Gen selection" section
    # ##########################################################
    # rdf_g.add_hists_singlecollection(dfGenP, histograms, sufElPt5)
    # add_puppicands_by_pdg(dfGenP, histograms, '')
    # ##########################################################

    # SPLIT ETA REGIONS BASED ON TKEL
    dfGenP = rdf_g.define_newcollection(dfGenP, sufElPt5, f'abs({sufElPt5}_eta) <= 1.47', 'EB')
    dfGenP = rdf_g.define_newcollection(dfGenP, sufElPt5, f'abs({sufElPt5}_eta) > 1.47 && abs({sufElPt5}_eta) <= 2.5', 'EE')

    gen_dRcuts = {'EB': 0.04, 'EE': 0.05}

    for ERegion in ['EB', 'EE']:
    # for ERegion in ['EB']:
        sufGenDYP = f'{sufGen}_DYP'
        sufElPt5ER = f'{sufElPt5}_{ERegion}'
        dfGenER = dfGenP.Filter(f'{sufElPt5ER}_n > 0', f'DYP{ERegion}')
        ut.create_rdf_checkpint(dfGenP, dfGenER, f"Applying selection: > 0 TkEl in region {ERegion}...")

        # # STEP 2_0_0 AND 2_1_0: GEN MATCH BLOCK - PRE GEN MATCH
        # #########################################################
        # dfGenER = anut.angdiff_hists(dfGenER, sufGenDYP, sufElPt5ER)
        # rdf_g.add_hists_multiplecolls(dfGenER, histograms, [f'{sufElPt5ER}', f'{sufGenDYP}_{sufElPt5ER}'])
        # anut.add_genmatching_efficiency_with_dRcut(histograms, f'DYP{ERegion}_{sufGenDYP}_{sufElPt5ER}')
        # ##########################################################

        # GEN MATCH
        dfGenER = anut.do_gen_match(dfGenER, sufGenDYP, sufElPt5ER, gen_dRcuts[ERegion])
        dfGenER = rdf_g.define_newcollection(dfGenER, sufGenDYP, f'{sufGenDYP}_recoidx != -1', 'MCH')
        dfGenER = rdf_g.define_newcollection(dfGenER, sufElPt5ER, f'{sufElPt5ER}_genidx != -1', 'MCH')
        
        sufGenMch = f'{sufGenDYP}_MCH'
        sufElPt5Mch = f'{sufElPt5ER}_MCH'

        # # STEP 2_0_0: GEN MATCH BLOCK - POST GEN MATCH
        # #########################################################
        # dfGenER = anut.angdiff_hists(dfGenER, sufGenMch, sufElPt5Mch)
        # rdf_g.add_hists_singlecollection(dfGenER, histograms, f'{sufGenMch}_{sufElPt5Mch}')
        # ##########################################################

        # STEP 2_1_0 AND 3_0_0: For comparing all TkEl to gen-matched TkEl
        #########################################################
        rdf_g.add_hists_singlecollection(dfGenER, histograms, sufElPt5Mch)
        ana_man.add_dataframe(key=f'DYP{ERegion}', df=dfGenER)
        #########################################################

        # Filter for atleast one gen-match TkEl in the defined eta region
        dfGenMER = dfGenER.Filter(f'{sufElPt5Mch}_n > 0', f'GM')
        ut.create_rdf_checkpint(dfGenER, dfGenMER, f"Applying selection: > 0 Gen-matched TkEl in region {ERegion}")

        # STEP 3_0_0: Add charged contribution to iso calc for gen-matched TkEl
        #########################################################
        dfGenMER = reiso.recalculate_puppi_iso(dfGenMER, sufElPt5Mch, sufPu)
        # dfGenMER = reiso.recalculate_puppi_iso(dfGenMER, sufElPt5Mch, sufPuPt1, drminlist=[0.01], drmax=0.4, ptmin=2, dzmax=1.0)
        # dfGenMER.Describe().Print()
        ana_man.add_dataframe(key=f'DYPM{ERegion}', df=dfGenMER)
        rdf_g.add_hists_multiplecolls(dfGenMER, histograms, [sufElPt5Mch,
                                        sufElPt5Mch+r'_reisotot2026:dRmin\d_\d{1,2}',
                                        sufElPt5Mch+r'_reisotot:dRmin\d_\d{1,2}',
                                        sufElPt5Mch+r'_reisooth:dRmin\d_\d{1,2}',
                                        sufElPt5Mch+r'_reisochg:dRmin\d_\d{1,2}',
                                        sufElPt5Mch+r'_reisonut:dRmin\d_\d{1,2}'])
        #########################################################

        dfGenMER = anut.make_puppi_by_angdiff_from_tkel(dfGenER, sufElPt5Mch, histograms, refPu=sufPuPt1)

    ana_man.add_histograms(histograms)
    return ana_man

        # STEP 4_0_0: Compare in pT separate bins
        #########################################################
        # dfGenMER = rdf_g.define_newcollection(dfGenMER, sufElPt5Mch, f"{sufElPt5Mch}_pt > 50", "Pt50")
        # dfGenMER = rdf_g.define_newcollection(dfGenMER, sufElPt5Mch, f"{sufElPt5Mch}_pt > 20 && {sufElPt5Mch}_pt <= 50", "20Pt50")
        # rdf_g.add_hists_multiplecolls(dfGenMER, histograms, [f"{sufElPt5Mch}_Pt50", f"{sufElPt5Mch}_20Pt50"])
        #########################################################

        # dfGenMER.Describe().Print()

    # # # Checked for pT sorting of the TkEl collection
    # dfGenER = dfGenER.Define(f'{sufElPt5Mch}_pt_sorted', f'checksorting<float>({sufElPt5Mch}_pt, true)')
    # histograms.append( dfGenER.Histo1D((f'{sufElPt5Mch}_pt_sorted', 'pt_sort', 6, -2, 4), f'{sufElPt5Mch}_pt_sorted') )
    # cp_pt_sort = anut.check_histogram_for_value(histograms, f'{sufElPt5Mch}_pt_sorted', bin=3) # bin 3 = value 0
    # if cp_pt_sort:
    #     raise RuntimeError("The Track Electrons pT's are not sorted from highest to lowest")

    # For each TkEl, observe the feature of PuppiCandidates around the TkEl in annular dR window
    # add_puppicands_by_pdg(dfGenER, histograms, '', tkelobj=sufElPt5Mch)
    # dfGenER = anut.make_puppi_by_angdiff_from_tkel(dfGenER, sufElPt5Mch, histograms)

    # Separate the lead and sub-lead gen matched TkEl
    dfgEBel1 = dfgenEB.Filter(f'{sufElPt5}_MCH_n >= 1', 'genDYEBTkEl1')
    dfgEBel1 = rdf_g.define_newcollection(dfgEBel1, f'{sufElPt5}_MCH', '0', 'El1')
    rdf_g.add_hists_singlecollection(dfgEBel1, histograms, f'{sufElPt5}_MCH_El1')
    # Puppi isolation and components by pdgId
    dfgEBel1 = reiso.recalcpuppiiso_comps_oneel(dfgEBel1, f'{sufElPt5}_MCH_El1', sufPuPt1)
    rdf_g.add_hists_singlecollection(dfgEBel1, histograms, f'{sufElPt5}_MCH_El1_Re',
                                     'dRmin\\d_\\d{1,2}_dRmax\\d_\\d{1,2}_[a-z0-9]+')

    # dfgEBel2 = dfgenEB.Filter(f'{sufElPt5}_MCH_n >= 2', 'genDYEBTkEl2')
    # dfgEBel2 = rdf_g.define_newcollection(dfgEBel2, f'{sufElPt5}_MCH', '1', 'El2')
    # rdf_g.add_hists_singlecollection(dfgEBel2, histograms, f'{sufElPt5}_MCH_El2')

    # Calculate the invariant mass
    # dfgEBel2 = dfgEBel2.Define(f'{sufElPt5}_MCH_m', f"ROOT::VecOps::RVec<float>({sufElPt5}_MCH_pt.size(), 0.000511)")
    # dfgEBel2 = dfgEBel2.Define(f'{sufElPt5}_MCH_Mee', f'ROOT::VecOps::InvariantMass<float>({sufElPt5}_MCH_pt,\
    # {sufElPt5}_MCH_eta, {sufElPt5}_MCH_phi, {sufElPt5}_MCH_m)')
    # histograms.append( dfgEBel2.Histo1D((f'genDYEBTkEl2_{sufElPt5}_MCH_Mee', 'mass [GeV]', 201, -1, 200),\
    # f'{sufElPt5}_MCH_Mee') )

    # dfgenEB_tkel2.Describe().Print()

    return histograms
