import time
import warnings

from ROOT import RDataFrame

from rdf_generic import add_hists_multiplecolls
from varmetadata import linkvartohist as v2h

# Define suffixes for different collections
sufEl = 'TkEleL2'
sufGen = 'GenEl'
sufPu = 'L1PuppiCands'


# Decorator to measure the execution time of a function
def time_eval(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"Execution time of function {func.__name__}: {execution_time:.6f} seconds")
        return result
    return wrapper


# Order of the variables matters here
# Identical order as getminangs() in define_cpp_utils.py
angdiff_vars = ['deta', 'dphi', 'dR']


def angdiff_hists(df: RDataFrame, referencecoll: str, targetcoll: str):

    collectionkey = f'{referencecoll}_{targetcoll}'
    get_angdiffs_str = f'getminangs({referencecoll}_eta, {referencecoll}_phi,\
                                    {targetcoll}_eta, {targetcoll}_phi)'

    for i, var in enumerate(angdiff_vars):
        df = df.Define(f'{collectionkey}_{var}', f'std::get<{i}>({get_angdiffs_str})')

    return df


def do_gen_match(df: RDataFrame, gencoll: str, recocoll: str, dRcut: float):

    perform_genmatch_str = f'getmatchedidxs({gencoll}_eta, {gencoll}_phi,\
                                            {recocoll}_eta, {recocoll}_phi, {dRcut})'
    df = df.Define(f'{gencoll}_recoidx', f'std::get<0>({perform_genmatch_str})')
    df = df.Define(f'{recocoll}_genidx', f'std::get<1>({perform_genmatch_str})')

    return df


def add_puppicands_by_pdg(df, histograms, collsuf):
    collection = sufPu if collsuf == '' else f'{sufPu}_{collsuf}'

    dfcolumnnames = list(df.GetColumnNames())
    collkey_vars = [columnname.split('_')[-1] for columnname in dfcolumnnames if columnname.startswith(f'{collection}')]

    pdict = {'El': 11, 'Mu': 13, 'Ph': 22, 'Nh': 130, 'Ch': 211}
    for k, v in pdict.items():
        sel = f'abs({collection}_pdgId) == {v}'
        psuf = f'{sufPu}{k}' if collsuf == '' else f'{sufPu}{k}_{collsuf}'

        for var in collkey_vars:
            if ('RVec' in df.GetColumnType(f'{collection}_{var}')):
                df = df.Define(f'{psuf}_{var}', f'{collection}_{var}[{sel}]')
                # if var in v2h:
                #     histograms.append(df.Histo1D((f'{psuf}_{var}', v2h[var][3],
                #                 v2h[var][0], v2h[var][1], v2h[var][2]), f'{psuf}_{var}'))
        if df.GetColumnType(f'{psuf}_pt').startswith('ROOT::VecOps::RVec'):
            df = df.Define(f'{psuf}_n', f'{psuf}_pt.size()')
            # histograms.append(df.Histo1D((f'{psuf}_n', v2h['n'][3],
            #                     v2h['n'][0], v2h['n'][1], v2h['n'][2]), f'{psuf}_n'))

    add_hists_multiplecolls(df, histograms, [f'{sufPu}{k}:{collsuf}' for k in pdict.keys()])


def add_genmatching_efficiency_with_dRcut(histograms, coll):
    orighistobj = [h for h in histograms if h.GetName() == f'{coll}_dR']
    hist = orighistobj[0].Clone()

    uflow = hist.GetBinContent(0)
    oflow = hist.GetBinContent(hist.GetNbinsX()+1)
    integral = hist.Integral()+uflow+oflow
    if integral != 0:
        hist.Scale(1.0 / integral)
        hist_integral = hist.GetCumulative(True, 'cumulative')
        binc = [hist_integral.GetBinContent(i) for i in range(hist_integral.GetNbinsX()+1)]
        bine = [hist_integral.GetBinError(i) for i in range(hist_integral.GetNbinsX()+1)]
        nseek = int(0.1*hist_integral.GetNbinsX())
        bindec = [binc[i]+bine[i] > binc[i+nseek] for i in range(hist_integral.GetNbinsX()+1-nseek)]
        stablebinpos = next((i for i, v in enumerate(bindec) if v), -1)
        print(f'In hist {hist.GetName()}, efficiency stabilises for:'\
              f' {hist_integral.GetBinCenter(stablebinpos):.4f} at {hist_integral.GetBinContent(stablebinpos):.4f}')
        histograms.append(hist_integral)
    else:
        warnings.warn(f"Integral for hist {hist.GetName()} null. Unable to produce gen matching efficiency with dR cut")


def check_histogram_for_value(histograms, histname, *, bin):
    hist = [h for h in histograms if h.GetName() == histname]
    if len(hist) > 0:
        hist = hist[0]
    if(hist.GetBinContent(bin) != 0):
        return True

    return False