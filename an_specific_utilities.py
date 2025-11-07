import array
import math as m
import numpy as np
import time
import warnings

from ROOT import RDataFrame

import rdf_generic as rdf_g
from rdf_generic import add_hists_multiplecolls, define_newcollection
from varmetadata import linkvartohist as v2h

# Define suffixes for different collections
sufEl = 'TkEleL2'
sufGen = 'GenEl'
sufPu = 'L1PuppiCands'


customisation_conds = {
    'canvas': [
        (lambda histname: 'Iso' in histname, lambda canvas: canvas.SetLogx()),
        # (lambda histname: 'bin2dR' in histname, lambda canvas: canvas.SetLogx())
    ]
}


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


def add_puppicands_by_pdg(df, histograms, collsuf, *, tkelobj=None):
    collection = sufPu if collsuf == '' else f'{sufPu}_{collsuf}'

    # if tkelobj is not None:
    #     get_angdiffs_str = f'getmindR_PuppiCandSize_TkElRef({collection}_eta, {collection}_phi, {collection}_pdgId,\
    #                                     {tkelobj}_eta, {tkelobj}_phi, {tkelobj}_caloEta, {tkelobj}_caloPhi)'
    #     df = df.Define(f'{collection}_bin2dR', get_angdiffs_str)

    dfcolumnnames = list(df.GetColumnNames())
    collkey_vars = [columnname.split('_')[-1] for columnname in dfcolumnnames if columnname.startswith(f'{collection}')]

    pdict = {'El': 11, 'Mu': 13, 'Ph': 22, 'Nh': 130, 'Ch': 211}
    for k, v in pdict.items():
        sel = f'abs({collection}_pdgId) == {v}'
        psuf = f'{sufPu}{k}' if collsuf == '' else f'{sufPu}{k}_{collsuf}'

        for var in collkey_vars:
            if ('RVec' in df.GetColumnType(f'{collection}_{var}')):
                df = df.Define(f'{psuf}_{var}', f'{collection}_{var}[{sel}]')

        if df.GetColumnType(f'{psuf}_pt').startswith('ROOT::VecOps::RVec'):
            df = df.Define(f'{psuf}_n', f'{psuf}_pt.size()')

        if tkelobj is not None:
            get_angdiffs_str = f"getminangs({psuf}_eta, {psuf}_phi, {tkelobj}_eta, {tkelobj}_phi)"
            if v == 22 or v == 130:
                get_angdiffs_str = f"getminangs({psuf}_eta, {psuf}_phi, {tkelobj}_caloEta, {tkelobj}_caloPhi)"
            df = df.Define(f'{psuf}_bin2dR', f"std::get<2>({get_angdiffs_str})")

    rdf_g.add_hists_multiplecolls(df, histograms, [f'{sufPu}{k}:{collsuf}' for k in pdict.keys()])

    if tkelobj is not None:
        for k, v in pdict.items():
            psuf = f'{sufPu}{k}' if collsuf == '' else f'{sufPu}{k}_{collsuf}'
            dfcolumnnames = list(df.GetColumnNames())
            collkey_vars = [columnname.split('_')[-1] for columnname in dfcolumnnames if columnname.startswith(psuf)]

            dRcond = f"{psuf}_bin2dR < 0.5"
            psuf_tkel = f'{sufPu}InTkEl{k}' if collsuf == '' else f'{sufPu}InTkEl{k}_{collsuf}'
            for var in collkey_vars:
                if ('RVec' in df.GetColumnType(f'{psuf}_{var}')):
                    df = df.Define(f'{psuf_tkel}_{var}', f'{psuf}_{var}[{dRcond}]')
            if df.GetColumnType(f'{psuf_tkel}_pt').startswith('ROOT::VecOps::RVec'):
                df = df.Define(f'{psuf_tkel}_n', f'{psuf_tkel}_pt.size()')

            duplication_dRcond = f"{psuf}_bin2dR > 0.005 && {psuf}_bin2dR < 0.5"
            psuf_nodup = f'{sufPu}NoTkElD{k}' if collsuf == '' else f'{sufPu}NoTkElD{k}_{collsuf}'
            for var in collkey_vars:
                if ('RVec' in df.GetColumnType(f'{psuf}_{var}')):
                    df = df.Define(f'{psuf_nodup}_{var}', f'{psuf}_{var}[{duplication_dRcond}]')
            if df.GetColumnType(f'{psuf_nodup}_pt').startswith('ROOT::VecOps::RVec'):
                df = df.Define(f'{psuf_nodup}_n', f'{psuf_nodup}_pt.size()')

        rdf_g.add_hists_multiplecolls(df, histograms, [f'{sufPu}InTkEl{k}:{collsuf}' for k in pdict.keys()])
        rdf_g.add_hists_multiplecolls(df, histograms, [f'{sufPu}NoTkElD{k}:{collsuf}' for k in pdict.keys()])


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


ndRbins = 10
ang_lsb = m.pi/720
d_ang = np.arange(ndRbins)*ang_lsb
fw_dRbins = [round(m.sqrt(e + p), 6) for e in d_ang**2 for p in d_ang**2]
fw_dRbins = sorted(list(set(fw_dRbins)))


def make_puppi_by_angdiff_from_tkel(df, refEg, histograms):
    dRs = [[0.0, 0.1], [0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 10]]
    # dRs = [[0.2, 0.3]]

    for dRrange in dRs:
        dRmin = dRrange[0]
        dRmax = dRrange[1]
        dRstr = f'{str(dRmin).replace('.', 'p')}dR{str(dRmax).replace('.', 'p')}'

        getpuppimask_annulardr_str = f"getpuppimask_annulardR({sufPu}_eta, {sufPu}_phi, {sufPu}_pdgId,\
            {refEg}_eta, {refEg}_caloEta, {refEg}_phi, {refEg}_caloPhi, {dRmin}, {dRmax})"
    
        df = df.Define(f'{sufPu}_mask{dRstr}', getpuppimask_annulardr_str)

        df = rdf_g.define_newcollection(df, sufPu, f'{sufPu}_mask{dRstr} == 1', dRstr)
        add_puppicands_by_pdg(df, histograms, dRstr)

    return df


def conditionally_modify_plots(histlist):

    newhistlist = []
    
    for hist in histlist:
        if hist.GetName().endswith('bin2dR'):
            bin2dR_edges = np.r_[0, np.arange(1, 10)*1e-4, np.arange(1, 10)*1e-3,
                                 np.arange(10, 100, 2)*1e-3, 0.1, 0.14, 0.2, 0.3, 0.5].tolist()

            hist = hist.Rebin(len(bin2dR_edges)-1, hist.GetName(), array.array('d', bin2dR_edges))
            hist.Scale(1.0, "width")
        
        newhistlist.append(hist)

    return newhistlist