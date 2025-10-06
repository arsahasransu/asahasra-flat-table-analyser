from ROOT import RDataFrame

import calc_puppi_iso as reiso
import rdf_generic as rdf_g
import an_specific_utilities as ut
from an_specific_utilities import sufEl, sufPu
from an_specific_utilities import add_puppicands_by_pdg


@ut.time_eval
def qcd_ana_main(df: RDataFrame):

    histograms = []

    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    add_puppicands_by_pdg(df, histograms, '')
    # df = ut.make_puppi_by_angdiff_from_tkel(df, sufEl, histograms)
    
    df = rdf_g.define_newcollection(df, sufEl, f'abs({sufEl}_eta) <= 1.4', 'EB')
    df = rdf_g.define_newcollection(df, sufEl, f'abs({sufEl}_eta) > 1.4 && abs({sufEl}_eta) <= 1.6', 'EM')
    df = rdf_g.define_newcollection(df, sufEl, f'abs({sufEl}_eta) > 1.6 && abs({sufEl}_eta) <= 2.1', 'EE')
    df = rdf_g.define_newcollection(df, sufEl, f'abs({sufEl}_eta) > 2.1', 'EF')

    for ERegion in ['EB', 'EM', 'EE', 'EF']:
    # for ERegion in ['EB']:
        sufElER = f'{sufEl}_{ERegion}'
        dfER = df.Filter(f'{sufElER}_n > 0', f'el{ERegion}')

        # rdf_g.add_hists_multiplecolls(df, histograms, [sufElER, sufPu])
        add_puppicands_by_pdg(dfER, histograms, '', tkelobj=sufElER)
        rdf_g.add_hists_multiplecolls(dfER, histograms, [sufElER])

        # dfER = ut.make_puppi_by_angdiff_from_tkel(dfER, sufElER, histograms)

    return histograms

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

