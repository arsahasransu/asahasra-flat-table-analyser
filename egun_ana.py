from ROOT import RDataFrame

import an_specific_utilities as anut
from an_specific_utilities import sufEl, sufGen, sufPu
from an_specific_utilities import add_puppicands_by_pdg
import pypkg.calc_puppi_iso as reiso
import pypkg.my_py_generic_utils as ut
import rdf_generic as rdf_g


@ut.time_eval
def egun_ana_main(ana_man: anut.SampleRDFManager) -> anut.SampleRDFManager:

    histograms = []

    df = ana_man.parent_df
    df = df.Define(sufGen+'_n', sufGen+'_pt.size()')
    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    df = df.Define(sufEl+'_absTkIso', f'{sufEl}_pt*{sufEl}_tkIso')

    # STEP 1_0_0: Enable for plots in "Gen properties" section
    ##########################################################
    # Confirmed that all gen electrons have status prompt == 2
    rdf_g.add_hists_multiplecolls(df, histograms, [sufGen, sufEl, sufPu])
    ##########################################################

    df = rdf_g.define_newcollection(df, sufGen, f'{sufGen}_prompt==2 && abs({sufGen}_eta)<=2.5', 'DYP')
    df = rdf_g.define_newcollection(df, sufEl, f'{sufEl}_pt >= 5.0', 'Pt5')
    sufElPt5 = sufEl+'_Pt5'
    df = rdf_g.define_newcollection(df, sufPu, f'{sufPu}_pt >= 1.0', 'Pt1')
    sufPuPt1 = sufPu+'_Pt1'
    rdf_g.add_hists_singlecollection(df, histograms, f'{sufGen}_DYP')

    dfGenP = df.Filter(f'{sufGen}_DYP_n > 0 && {sufElPt5}_n > 0', 'genDYP')
    ana_man.add_dataframe(key='genDYP', df=dfGenP)
    ut.create_rdf_checkpint(df, dfGenP, "Applying selection: >0 Gen electron with prompt "
                                        "status 2 and |eta| < 2.5... \n>0 reconstructed TkEl...")

    # SPLIT ETA REGIONS BASED ON TKEL
    dfGenP = rdf_g.define_newcollection(dfGenP, sufElPt5, f'abs({sufElPt5}_eta) <= 1.47', 'EB')
    dfGenP = rdf_g.define_newcollection(dfGenP, sufElPt5, f'abs({sufElPt5}_eta) > 1.47 && abs({sufElPt5}_eta) <= 2.5', 'EE')

    gen_dRcuts = {'EB': 0.03, 'EE': 0.04}

    for ERegion in ['EB', 'EE']:
    # for ERegion in ['EB']:
        sufGenDYP = f'{sufGen}_DYP'
        sufElPt5ER = f'{sufElPt5}_{ERegion}'
        dfGenER = dfGenP.Filter(f'{sufElPt5ER}_n > 0', f'DYP{ERegion}')
        ut.create_rdf_checkpint(dfGenP, dfGenER, f"Applying selection: > 0 TkEl in region {ERegion}...")

        # STEP 2_0_0 AND 2_1_0: GEN MATCH BLOCK - PRE GEN MATCH
        #########################################################
        dfGenER = anut.angdiff_hists(dfGenER, sufGenDYP, sufElPt5ER)
        rdf_g.add_hists_multiplecolls(dfGenER, histograms, [f'{sufElPt5ER}', f'{sufGenDYP}_{sufElPt5ER}'])
        anut.add_genmatching_efficiency_with_dRcut(histograms, f'DYP{ERegion}_{sufGenDYP}_{sufElPt5ER}')
        ##########################################################

        # GEN MATCH
        dfGenER = anut.do_gen_match(dfGenER, sufGenDYP, sufElPt5ER, gen_dRcuts[ERegion])
        dfGenER = rdf_g.define_newcollection(dfGenER, sufGenDYP, f'{sufGenDYP}_recoidx != -1', 'MCH')
        dfGenER = rdf_g.define_newcollection(dfGenER, sufElPt5ER, f'{sufElPt5ER}_genidx != -1', 'MCH')
        
        sufGenMch = f'{sufGenDYP}_MCH'
        sufElPt5Mch = f'{sufElPt5ER}_MCH'

        # STEP 2_0_0: GEN MATCH BLOCK - POST GEN MATCH
        #########################################################
        dfGenER = anut.angdiff_hists(dfGenER, sufGenMch, sufElPt5Mch)
        rdf_g.add_hists_singlecollection(dfGenER, histograms, f'{sufGenMch}_{sufElPt5Mch}')
        ##########################################################

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

        dfGenMERP = anut.make_puppi_by_angdiff_from_tkel(dfGenMER, sufElPt5Mch, histograms, refPu=sufPuPt1)
        ana_man.add_dataframe(key=f'DYPM{ERegion}P', df=dfGenMERP)

    ana_man.add_histograms(histograms)
    return ana_man
