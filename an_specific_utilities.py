import array
import math as m
import numpy as np
import warnings

from ROOT import RDataFrame, TFile

import my_py_generic_utils as ut
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
        # (lambda histname: histname.endswith('eta') or histname.endswith('Eta'), lambda canvas: canvas.SetLogy(False))
        # (lambda histname: 'bin2dR' in histname, lambda canvas: canvas.SetLogx())
    ]
}


class SampleRDFManager:
    def __init__(self, s_name:str, s_opts:dict, *,
                 opts = None,
                 tree_name: str = "Events"):
        self.s_name = s_name
        self.file_pattern = (opts['input_dir_prefix'] if opts is not None else '')+'/'+s_opts['input_file_pattern']
        self._opts = opts
        if self._opts is not None:
            self.tree_name = self._opts.get('tree_name', tree_name)
        else:
            self.tree_name = tree_name
        self.max_events = s_opts.get('max_events', 0)

        self.parent_df = self.load_rdf()
        self.df_dict = {}
        self.histograms = []

    def load_rdf(self):
        try:
            df = RDataFrame(self.tree_name, self.file_pattern)
            print("\n")
            print(f"PROCESSING {df.Count().GetValue()} EVENTS IN SAMPLE {self.s_name.upper()}")
            return df
        except FileNotFoundError as err:
            raise FileNotFoundError(f"Could not find files for sample {self.s_name} in path {self.file_pattern}: {err}")
        except Exception as e:
            raise Exception(f"An error occurred while loading the RDF for sample {self.s_name}: {e}")
    
    def add_dataframe(self, key: str, df: RDataFrame):
        self.df_dict[key] = df

    def get_dataframe(self, key: str) -> RDataFrame:
        if key not in self.df_dict:
            raise KeyError(f"DataFrame with key '{key}' not found.")
        
        return self.df_dict.get(key)
    
    def add_histograms(self, histolist: list):
        self.histograms.extend(histolist)

    def get_histograms(self) -> list:
        return self.histograms


# Identical order as getminangs() in define_cpp_utils.py
angdiff_vars = ['deta', 'dphi', 'dR']


def angdiff_hists(df: RDataFrame, referencecoll: str, targetcoll: str):

    collectionkey = f'{referencecoll}_{targetcoll}'
    get_angdiffs_str = f'getminangs({referencecoll}_eta, {referencecoll}_phi,\
                                    {targetcoll}_eta, {targetcoll}_phi)'

    for i, var in enumerate(angdiff_vars):
        df = df.Define(f'{collectionkey}_{var}', f'{get_angdiffs_str}.at("{var}")')

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
        nseek = int(0.05*hist_integral.GetNbinsX())
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


def make_plotnorm_by_scheme(hcoll_2d, scheme, *, summed_sample_pos=None, byevent_vartag='n'):

    if scheme == "default_no_norm":
        return [[1.0 for h in hcoll_1d] for hcoll_1d in hcoll_2d]

    if scheme == "hist_integral":
        return [[h.Integral(0, h.GetNbinsX()+1) for h in hcoll_1d] for hcoll_1d in hcoll_2d]

    if scheme == "summed_components":
        if summed_sample_pos == None:
            raise ValueError("Reference summed sample index not provided to normalise with \"summed_componenets\" strategy!")
        if summed_sample_pos > len(hcoll_2d[0]):
            raise ValueError("Reference summed sample index larger than the length of the histogram list for each variable."\
                             "Potentially forgot to deduct 1 for array indexing or "\
                             "accumulated histograms along the wrong dimension of the histogram list.")
        ssp = summed_sample_pos
        f_nh = lambda h: h.Integral(0, h.GetNbinsX()+1)

        check_byevent_vartag_collection = False
        notag_nev = []
        nev = []
        for hcoll_1d in hcoll_2d:
            name_h1d = [h.GetName() for h in hcoll_1d]
            compare_vartag = [n.endswith("_"+byevent_vartag) for n in name_h1d]

            if all(compare_vartag):
                notag_nev = [f_nh(h) for i,h in enumerate(hcoll_1d) if i != ssp]
                nev = [(f_nh(h) if i != ssp else -1.0) for i,h in enumerate(hcoll_1d)]
                check_byevent_vartag_collection = True
                break

        if not check_byevent_vartag_collection:
            raise AttributeError("Did not find the variable for event count while "\
                                 "trying \"summed_component\" normalisation.")

        notag_nh = [[f_nh(h) for i,h in enumerate(h1d) if i != ssp] for h1d in hcoll_2d]
        notag_normh = [[nh/ne for nh,ne in zip(nh1var, notag_nev)] for nh1var in notag_nh]
        notag_norms = [sum(normh1var) for normh1var in notag_normh]

        norm = []
        for notag_norm,hcoll_1d in zip(notag_norms,hcoll_2d):
            name_h1d = [h.GetName() for h in hcoll_1d]
            compare_vartag = [name.endswith("_"+byevent_vartag) for name in name_h1d]

            if not all(compare_vartag):
                norm.append([(nev[i]*notag_norm if i != ssp else f_nh(h)) for i,h in enumerate(hcoll_1d)])
            else:
                norm.append([f_nh(h) for h in hcoll_1d])

        return norm


@ut.time_eval
def writehists_to_file(ofname, histograms):
    outfile = TFile(ofname, 'RECREATE')
    for hist in histograms:
        hist.Write()
    outfile.Close()
