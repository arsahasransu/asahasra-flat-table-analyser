from ROOT import RDataFrame
import ROOT

import an_specific_utilities as anut
from an_specific_utilities import sufEl, sufPu
from an_specific_utilities import add_puppicands_by_pdg
import calc_puppi_iso as reiso
import my_py_generic_utils as ut
import rdf_generic as rdf_g


@ut.time_eval
def qcd_ana_main(ana_man: anut.SampleRDFManager) -> anut.SampleRDFManager:

    histograms = []

    df = ana_man.parent_df
    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    df = rdf_g.define_newcollection(df, sufEl, f'{sufEl}_pt > 5', 'Pt5')
    sufElPt5 = sufEl+'_Pt5'

    df = df.Define(sufElPt5+'_absTkIso', f'{sufElPt5}_pt*{sufElPt5}_tkIso')
    
    dfE = df.Filter(f'{sufElPt5}_n > 0', 'tke1')
    ut.create_rdf_checkpint(df, dfE, "Applying selection: >0 reconstructed TkEl...")

    # STEP 1_2_0: Enable for plots in "Gen selection" section
    ##########################################################
    # rdf_g.add_hists_singlecollection(dfE, histograms, sufElPt5)
    # add_puppicands_by_pdg(dfE, histograms, '')
    ##########################################################
    # df = anut.make_puppi_by_angdiff_from_tkel(df, sufElPt5, histograms)
    
    dfE = rdf_g.define_newcollection(dfE, sufElPt5, f'abs({sufElPt5}_eta) <= 1.5', 'EB')
    dfE = rdf_g.define_newcollection(dfE, sufElPt5, f'abs({sufElPt5}_eta) > 1.5 && abs({sufElPt5}_eta) <= 2.5', 'EE')

    for ERegion in ['EB', 'EE']:
    # for ERegion in ['EB']:
        sufElPt5ER = f'{sufElPt5}_{ERegion}'
        dfER = dfE.Filter(f'{sufElPt5ER}_n > 0', ERegion)
        ut.create_rdf_checkpint(dfE, dfER, f"Applying selection: > 0 TkEl in region {ERegion}...")

        # STEP 3_0_0: sufElER necessary to make reference ROC curves for tkIso
        #########################################################
        dfER = reiso.recalculate_puppi_iso(dfER, sufElPt5ER, sufPu)
        ana_man.add_dataframe(key=ERegion, df=dfER)
        rdf_g.add_hists_multiplecolls(dfER, histograms, [sufElPt5ER,
                                        sufElPt5ER+'_reisotot:dRmin\\d_\\d{1,2}_[a-z0-9]+'])
        #########################################################

        # STEP 4_0_0: Compare in pT separate bins
        #########################################################
        # dfER = rdf_g.define_newcollection(dfER, sufElPt5ER, f"{sufElPt5ER}_pt > 50", "Pt50")
        # dfER = rdf_g.define_newcollection(dfER, sufElPt5ER, f"{sufElPt5ER}_pt > 20 && {sufElPt5ER}_pt <= 50", "20Pt50")
        # rdf_g.add_hists_multiplecolls(dfER, histograms, [f"{sufElPt5ER}_Pt50", f"{sufElPt5ER}_20Pt50"])
        #########################################################
    ana_man.add_histograms(histograms)
    return ana_man
