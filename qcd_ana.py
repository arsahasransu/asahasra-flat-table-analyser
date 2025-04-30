import utilities as utilities
from utilities import puppi_vars, tkell2_vars


sufEl = 'TkEleL2'
sufPu = 'L1PuppiCands'

def add_tkel_filter(df, prefiltid='', filterid='', filter_cond=None):

    precoll = f'{sufEl}{prefiltid}'
    collection = f'{sufEl}{filterid}'
    for k, _ in tkell2_vars.items():
        df = df.Define(f'{collection}_{k}', f'{precoll}_{k}[{filter_cond}]')
    df = df.Define(f'{collection}_n', f'{collection}_pt.size()')

    df = df.Filter(f'{collection}_n > 0')

    precoll = f'{sufPu}{prefiltid}'
    collection = f'{sufPu}{filterid}'
    for k, _ in puppi_vars.items():
        df = df.Define(f'{collection}_{k}', f'{precoll}_{k}')
    df = df.Define(f'{collection}_n', f'{collection}_pt.size()')

    return df


def add_el_hists(df, histograms, selid='', filterid=''):

    collection = f'{sufEl}{selid}'
    for k, (bins, low, up) in tkell2_vars.items():
        histograms.append(df.Histo1D((f'{collection}{filterid}_{k}', f'{k}', bins, low, up), 
                                     f'{collection}_{k}'))
    histograms.append(df.Histo1D((f'{collection}{filterid}_n', 'mult', 21, -1, 20), f'{collection}_n'))


def add_puppi_hists(df, histograms, selid='', filterid=''):

    collection = f'{sufPu}{selid}'
    for k, (bins, low, up) in puppi_vars.items():
        histograms.append(df.Histo1D((f'{collection}{filterid}_{k}', f'{k}', bins, low, up), 
                                     f'{collection}_{k}'))
    histograms.append(df.Histo1D((f'{collection}{filterid}_n', 'mult', 21, -1, 20), f'{collection}_n'))


def add_standard_hists(df, histograms, selid='', filterid=''):

    add_el_hists(df, histograms, selid, filterid)
    add_puppi_hists(df, histograms, selid, filterid)


@utilities.time_eval
def qcd_ana_main(df):

    histograms = []
    
    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    # Consistency check for pt sorting
    # TODO: Implement this

    # No selection
    add_standard_hists(df, histograms, '', '')

    df_TkElEB = add_tkel_filter(df, '', 'TkElPt10EB', f'{sufEl}_pt > 10 && abs({sufEl}_eta) < 1.479')
    add_standard_hists(df_TkElEB, histograms, 'TkElPt10EB', '')

    return (df, histograms)