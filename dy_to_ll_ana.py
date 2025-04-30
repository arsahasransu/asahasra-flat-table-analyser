import utilities as utilities
from utilities import genel_vars, genmch_vars, puppi_vars, tkell2_vars


sufGen = 'GenEl'
sufEl = 'TkEleL2'
sufPu = 'L1PuppiCands'

def do_gen_match(df, histograms, gencoll, recocoll):

    perform_genmatch_str = f'getmatchedidxs({gencoll}_eta, {gencoll}_phi,\
                                            {recocoll}_eta, {recocoll}_phi, 0.02)'
    df = df.Define(f'{gencoll}_reco_idx', f'std::get<0>({perform_genmatch_str})')
    df = df.Define(f'{recocoll}_gen_idx', f'std::get<1>({perform_genmatch_str})')

    for k, (bins, low, up) in genel_vars.items():
        df = df.Define(f'{gencoll}_{recocoll}matched_{k}', f'{gencoll}_{k}[{gencoll}_reco_idx != -1]')
        histograms.append(df.Histo1D((f'{gencoll}_{recocoll}matched_{k}', f'{k}', bins, low, up), 
                                      f'{gencoll}_{recocoll}matched_{k}'))
    df = df.Define(f'{gencoll}_{recocoll}matched_n', f'{gencoll}_{recocoll}matched_pt.size()')
    histograms.append(df.Histo1D((f'{gencoll}_{recocoll}matched_n', 'mult', 21, -1, 20),
                                  f'{gencoll}_{recocoll}matched_n'))
        
    for k, (bins, low, up) in tkell2_vars.items():
        df = df.Define(f'{recocoll}_{gencoll}matched_{k}', f'{recocoll}_{k}[{recocoll}_gen_idx != -1]')
        histograms.append(df.Histo1D((f'{recocoll}_{gencoll}matched_{k}', f'{k}', bins, low, up), 
                                      f'{recocoll}_{gencoll}matched_{k}'))
    df = df.Define(f'{recocoll}_{gencoll}matched_n', f'{recocoll}_{gencoll}matched_pt.size()')
    histograms.append(df.Histo1D((f'{recocoll}_{gencoll}matched_n', 'mult', 21, -1, 20),
                                  f'{recocoll}_{gencoll}matched_n'))

    return df


def pre_gen_mch_hists(df, histograms, gencoll, recocoll):

    genmchcoll = f'{gencoll}_{recocoll}'
    get_angdiffs_str = f'getminangs({gencoll}_eta, {gencoll}_phi,\
                                    {recocoll}_eta, {recocoll}_phi)'

    for i, (k, (bins, low, up)) in enumerate(genmch_vars.items()):
        df = df.Define(f'{genmchcoll}_{k}', f'std::get<{i}>({get_angdiffs_str})')
        histograms.append(df.Histo1D((f'{genmchcoll}_{k}', f'{k}', bins, low, up), 
                                     f'{genmchcoll}_{k}'))

    return df


def filter_by_gen(df):

    sufGenFilt = 'genEB'

    for k, _ in genel_vars.items():
        df = df.Define(f'{sufGen}{sufGenFilt}_{k}', f'{sufGen}_{k}[abs(GenEl_eta)<1.47 && GenEl_prompt==2]')
    df = df.Define(f'{sufGen}{sufGenFilt}_n', f'{sufGen}{sufGenFilt}_pt.size()')

    df = df.Filter(f'{sufGen}{sufGenFilt}_pt.size() >= 1')

    for k, _ in tkell2_vars.items():
        df = df.Define(f'{sufEl}{sufGenFilt}_{k}', f'{sufEl}_{k}')
    df = df.Define(f'{sufEl}{sufGenFilt}_n', f'{sufEl}{sufGenFilt}_pt.size()')

    for k, _ in puppi_vars.items():
        df = df.Define(f'{sufPu}{sufGenFilt}_{k}', f'{sufPu}_{k}')
    df = df.Define(f'{sufPu}{sufGenFilt}_n', f'{sufPu}{sufGenFilt}_pt.size()')

    return df


def add_gen_hists(df, histograms, selid='', filterid=''):

    collection = f'{sufGen}{selid}'
    for k, (bins, low, up) in genel_vars.items():
        histograms.append(df.Histo1D((f'{collection}{filterid}_{k}', f'{k}', bins, low, up), 
                                     f'{collection}_{k}'))
    histograms.append(df.Histo1D((f'{collection}{filterid}_n', 'mult', 21, -1, 20), f'{collection}_n'))


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

    add_gen_hists(df, histograms, selid, filterid)
    add_el_hists(df, histograms, selid, filterid)
    add_puppi_hists(df, histograms, selid, filterid)


@utilities.time_eval
def dy_to_ll_ana_main(df):

    histograms = []
    
    df = df.Define(sufGen+'_n', sufGen+'_pt.size()')
    df = df.Define(sufEl+'_n', sufEl+'_pt.size()')
    df = df.Define(sufPu+'_n', sufPu+'_pt.size()')

    # Consistency check for pt sorting
    # TODO: Implement this

    # No selection
    add_standard_hists(df, histograms, '', '')

    # Filter for gen electrons in EB
    dfEB = filter_by_gen(df)
    add_standard_hists(dfEB, histograms, 'genEB', '')

    # Observe angles before gen matching
    dfEB = pre_gen_mch_hists(dfEB, histograms, sufGen+'genEB', sufEl)
    dfEB = do_gen_match(dfEB, histograms, sufGen+'genEB', sufEl)

    return (df, histograms)