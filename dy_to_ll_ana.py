import utilities as utilities
from utilities import genel_vars, puppiel_vars, genmch_vars


def gen_match(df, gencoll, pupcoll):

    histograms = []

    perform_genmatch_str = f'perform_genmatch({gencoll}_pt, {gencoll}_eta, {gencoll}_phi,\
                                              {pupcoll}_pt, {pupcoll}_eta, {pupcoll}_phi, 0.1)'
    df = df.Define(f'{gencoll}_puppi_idx', f'std::get<0>({perform_genmatch_str})')
    df = df.Define(f'{pupcoll}_gen_idx', f'std::get<1>({perform_genmatch_str})')

    for k, (bins, low, up) in genel_vars.items():
        df = df.Define(f'{gencoll}_{pupcoll}matched_{k}', f'{gencoll}_{k}[{gencoll}_puppi_idx != -1]')
        histograms.append(df.Histo1D((f'{gencoll}_{pupcoll}matched_{k}', f'{k}', bins, low, up), 
                                     f'{gencoll}_{pupcoll}matched_{k}'))

    return (df, histograms)


def pre_gen_mch_hists(df, gencoll, pupcoll):

    histograms = []

    genmchcoll = f'{gencoll}_{pupcoll}'
    get_angdiffs_str = f'calculate_angdiff({gencoll}_pt, {gencoll}_eta, {gencoll}_phi, {gencoll}_charge,\
                                           {pupcoll}_pt, {pupcoll}_eta, {pupcoll}_phi, {pupcoll}_charge)'
    ctr = 0

    for k, (bins, low, up) in genmch_vars.items():
        df = df.Define(f'{genmchcoll}_{k}', f'std::get<{ctr}>({get_angdiffs_str})')
        histograms.append(df.Histo1D((f'{genmchcoll}_{k}', f'{k}', bins, low, up), 
                                     f'{genmchcoll}_{k}'))
        ctr += 1

    return (df, histograms)


def filter_by_gen(df):

    for k, _ in genel_vars.items():
        df = df.Define(f'GenElPromptEB_{k}', f'GenEl_{k}[abs(GenEl_eta)<1.48 && GenEl_prompt==2]')

    df = df.Define('GenElPromptEB_n', 'GenElPromptEB_pt.size()')
    df = df.Filter('GenElPromptEB_pt.size() >= 1')

    return df


def add_gen_hists(df, selid='', filterid=''):

    histograms = []

    collection = f'GenEl{selid}'
    for k, (bins, low, up) in genel_vars.items():
        histograms.append(df.Histo1D((f'{collection}{filterid}_{k}', f'{k}', bins, low, up), 
                                     f'{collection}_{k}'))
    histograms.append(df.Histo1D((f'{collection}{filterid}_n', 'mult', 21, -1, 20), f'{collection}_n'))

    return (df,histograms)


def filter_by_puppiel(df, prevselid, selid, selstring):

    oldcoll = f'PuppiEl{prevselid}'
    newcoll = f'PuppiEl{selid}'
    for k, _ in puppiel_vars.items():
        df = df.Define(f'{newcoll}_{k}', f'{oldcoll}_{k}[{selstring}]')
    df = df.Define(f'{newcoll}_n', f'{newcoll}_pt.size()')

    return df


def add_el_hists(df, selid='', filterid=''):

    histograms = []

    collection = f'PuppiEl{selid}'
    for k, (bins, low, up) in puppiel_vars.items():
        histograms.append(df.Histo1D((f'{collection}{filterid}_{k}', f'{k}', bins, low, up), 
                                     f'{collection}_{k}'))
    histograms.append(df.Histo1D((f'{collection}{filterid}_n', 'mult', 21, -1, 20), f'{collection}_n'))
    return (df, histograms)


@utilities.time_eval
def dy_to_ll_ana_main(df):

    histograms = []

    df = df.Define('GenEl_n', 'GenEl_pt.size()')
    df = df.Define('PuppiEl_n', 'PuppiEl_pt.size()')

    # Consistency check for pt sorting
    # TODO: Implement this

    # No selection
    df, more_hists = add_gen_hists(df, '')
    histograms.extend(more_hists)
    df,more_hists = add_el_hists(df, '', '')
    histograms.extend(more_hists)

    # Filter for prompt gen electrons in EB
    df = filter_by_gen(df)
    df,more_hists = add_gen_hists(df, 'PromptEB', '')
    histograms.extend(more_hists)
    df,more_hists = add_el_hists(df, '', '_genPromptEB')
    histograms.extend(more_hists)

    # Add filter for Puppi electrons in EB
    df = filter_by_puppiel(df, '', 'EB', 'abs(PuppiEl_eta)<1.479')
    df = df.Filter('PuppiElEB_pt.size() >= 1')
    df,more_hists = add_gen_hists(df, 'PromptEB', '_puppiEB')
    histograms.extend(more_hists)
    df,more_hists = add_el_hists(df, 'EB', '_genPromptEB')
    histograms.extend(more_hists)

    # Observe angles before gen matching
    df,more_hists = pre_gen_mch_hists(df, 'GenElPromptEB', 'PuppiElEB')
    histograms.extend(more_hists)

    # Perform gen matching
    df,more_hists = gen_match(df, 'GenElPromptEB', 'PuppiElEB')
    histograms.extend(more_hists)

    return (df, histograms)