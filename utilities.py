import string
import time

from ROOT import RDataFrame

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


def angdiff_hists(df:RDataFrame, referencecoll:string, targetcoll:string):

    collectionkey = f'{referencecoll}_{targetcoll}'
    get_angdiffs_str = f'getminangs({referencecoll}_eta, {referencecoll}_phi,\
                                    {targetcoll}_eta, {targetcoll}_phi)'

    for i, var in enumerate(angdiff_vars):
        df = df.Define(f'{collectionkey}_{var}', f'std::get<{i}>({get_angdiffs_str})')

    return df


def do_gen_match(df:RDataFrame, gencoll:string, recocoll:string):

    perform_genmatch_str = f'getmatchedidxs({gencoll}_eta, {gencoll}_phi,\
                                            {recocoll}_eta, {recocoll}_phi, 0.02)'
    df = df.Define(f'{gencoll}_recoidx', f'std::get<0>({perform_genmatch_str})')
    df = df.Define(f'{recocoll}_genidx', f'std::get<1>({perform_genmatch_str})')

    return df


# threep_vars = ['pt', 'eta', 'phi']


# charged_threep_vars = threep_vars + ['charge', 'vz']


# genel_vars = charged_threep_vars + ['caloeta', 'calophi', 'prompt']


# tkell2_vars = charged_threep_vars + ['hwQual', 'tkPt', 'caloEta',
#                                      'tkEta', 'caloPhi', 'tkPhi',
#                                      'pfIso', 'puppiIso', 'tkIso']


# puppi_vars = threep_vars + ['mass', 'hwTkQuality', 'pdgId',
#                             'puppiWeight', 'z0']
