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

    dfE = df.Filter(f'{sufEl}_n > 0', 'tke1')
    ut.create_rdf_checkpint(df, dfE, "Applying selection: >0 reconstructed TkEl...")

    # STEP 1_2_0: Enable for plots in "Gen selection" section
    ##########################################################
    # rdf_g.add_hists_singlecollection(dfE, histograms, sufEl)
    # add_puppicands_by_pdg(dfE, histograms, '')
    ##########################################################
    # df = anut.make_puppi_by_angdiff_from_tkel(df, sufEl, histograms)
    
    dfE = rdf_g.define_newcollection(dfE, sufEl, f'abs({sufEl}_eta) <= 1.5', 'EB')
    dfE = rdf_g.define_newcollection(dfE, sufEl, f'abs({sufEl}_eta) > 1.5 && abs({sufEl}_eta) <= 2.5', 'EE')

    # for ERegion in ['EB', 'EE']:
    for ERegion in ['EB']:
        sufElER = f'{sufEl}_{ERegion}'
        dfER = dfE.Filter(f'{sufElER}_n > 0', ERegion)
        ut.create_rdf_checkpint(dfE, dfER, f"Applying selection: > 0 TkEl in region {ERegion}...")

        rdf_g.add_hists_singlecollection(dfER, histograms, sufElER)

        # rdf_g.add_hists_multiplecolls(df, histograms, [sufElER, sufPu])
        # add_puppicands_by_pdg(dfER, histograms, '', tkelobj=sufElER)
        # rdf_g.add_hists_multiplecolls(dfER, histograms, [sufElER])

        # dfER = anut.make_puppi_by_angdiff_from_tkel(dfER, sufElER, histograms)

    ana_man.add_histograms(histograms)
    return ana_man

    # df = rdf_g.define_newcollection(df, sufEl, f'{sufEl}_pt > 10 && abs({sufEl}_eta) < 1.479', 'TkElEBPt10')
    # dfTkElEBPt10 = df.Filter(f'{sufEl}_TkElEBPt10_n > 0', 'tkelEBpt10')
    # dfTkElEBPt10 = rdf_g.define_newcollection(dfTkElEBPt10, f'{sufEl}_TkElEBPt10', '0', 'El1')
    # rdf_g.add_hists_multiplecolls(dfTkElEBPt10, histograms, [sufPu,
    #                                                          f'{sufEl}_TkElEBPt10_El1',
    #                                                          f'{sufEl}_TkElEBPt10_El2'])

    # dfTkElEBPt10 = reiso.recalculate_puppi_iso(dfTkElEBPt10, f'{sufEl}_TkElEBPt10', sufPu)
    # rdf_g.add_hists_singlecollection(dfTkElEBPt10, histograms, f'{sufEl}_TkElEBPt10_Re', 'dRmin\\d_\\d{1,2}')

    # Puppi isolation and components by pdgId
    # dfTkElEBPt10 = reiso.recalcpuppiiso_comps_oneel(dfTkElEBPt10, f'{sufEl}_TkElEBPt10_El1', sufPu)
    # rdf_g.add_hists_singlecollection(dfTkElEBPt10, histograms, f'{sufEl}_TkElEBPt10_El1_Re',\
    #                                  'dRmin\\d_\\d{1,2}_dRmax\\d_\\d{1,2}_[a-z0-9]+')

