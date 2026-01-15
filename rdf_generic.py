import operator
import re

import awkward as ak
import numpy as np

from ROOT import RDataFrame
import ROOT

import my_py_generic_utils as ut
from varmetadata import linkvartohist


# Filter on a single collections with a list of potential variables
def define_newcollection(df: RDataFrame, collection: str, selection: str, sid: str):
    dfcolumnnames = list(df.GetColumnNames())
    collkey_vars = list(set([columnname.split('_')[-1] for columnname in dfcolumnnames if columnname.startswith(collection+'_')]))
    for var in collkey_vars:
        if (f'{collection}_{var}' in dfcolumnnames) and ('RVec' in df.GetColumnType(f'{collection}_{var}')):
            df = df.Define(f'{collection}_{sid}_{var}',
                           f'{collection}_{var}[{selection}]')

    if df.GetColumnType(f'{collection}_{sid}_pt').startswith('ROOT::VecOps::RVec'):
        df = df.Define(f'{collection}_{sid}_n', f'{collection}_{sid}_pt.size()')
    return df


# Add histograms for a collection of variables
# Input: optionally the selection ID of the collection
# Infers the filter ID and the list of available variables for the collection
def add_hists_singlecollection(df: RDataFrame, histograms: list, collection: str, sreg: str = ''):

    # Infer the filter ID of the dataframe
    fid = 'base_'
    # TODO LINE BELOW GIVES MORE TIME DELAY THAN EXPECTED - REASON UNKNOWN
    filternameslist = df.GetFilterNames()
    if len(filternameslist) != 0:
        fid = filternameslist[-1] + '_'

    # Infer the available variables for the collection and selection ID
    collkey = collection
    dfcolumnnames = list(df.GetColumnNames())

    if sreg != '':
        fullcolumnregex = f'({collection})_({sreg})_(\\w+)'
        selection_regex = re.compile(fullcolumnregex)
        dfcolumn_regmatchobj = [re.fullmatch(selection_regex, str(dfcolname)) for dfcolname in dfcolumnnames]
        histograms.extend([df.Histo1D((f'{fid}{m.group(0)}', linkvartohist[m.group(3)][3],
                           linkvartohist[m.group(3)][0], linkvartohist[m.group(3)][1], linkvartohist[m.group(3)][2]),
                           f'{m.group(0)}') for m in dfcolumn_regmatchobj if m is not None and m.group(3) in linkvartohist])
    else:
        collkey_vars = [columnname.split('_')[-1] for columnname in dfcolumnnames if columnname.startswith(collkey+'_')]
        histograms.extend([df.Histo1D((f'{fid}{collkey}_{var}', linkvartohist[var][3],
                                       linkvartohist[var][0], linkvartohist[var][1], linkvartohist[var][2]),
                                      f'{collkey}_{var}') for var in collkey_vars if var in linkvartohist])


def add_hists_multiplecolls(df: RDataFrame, histograms: list, collections: list):
    for collection in collections:
        if ':' in collection:
            collectionparts = collection.split(':')
            add_hists_singlecollection(df, histograms, collectionparts[0], collectionparts[1])
        else:
            add_hists_singlecollection(df, histograms, collection)


def save_rdf_snapshot_to_parquet(df: RDataFrame, cols: list[str], savename: str):
    print("Saving: " + ", ".join(cols) + f"to {savename}.parquet")

    np_dict = df.AsNumpy(cols)
    # Convert each column to Awkward and combine into a record array
    ak_fields = {col: ak.from_iter(np_dict[col]) for col in cols}
    ak_table = ak.zip(ak_fields, depth_limit=1)  # creates a table with named fields

    # Save as Parquet (efficient, columnar)
    ak.to_parquet(ak_table, f"{savename}.parquet")


@ut.time_eval
def roc_from_rdfs(signal_rdf: RDataFrame, background_rdf: RDataFrame, signal_ref_col: str,
                  background_ref_col:str,
                  signal_weight_col: str = None, backgound_weight_col: str = None, *,
                  ngrid: int = 20, title: str = "default", order = "higher"):
                  # Other vars to add the function above
                  # wt_col=None, refvarbins, ngrid

    # Expect that the reference column is a variable of a collection which can have multiple entries per row(event)
    # If at least one entry in the event passes then the event passes

    title = f'{signal_ref_col}_{background_ref_col}_roc' if title == "default" else title
    is_higher = order == "higher"

    s_cols = [signal_ref_col] + ([signal_weight_col] if signal_weight_col else [])
    b_cols = [background_ref_col] + ([backgound_weight_col] if backgound_weight_col else [])
    np_s = signal_rdf.AsNumpy(s_cols)
    np_b = background_rdf.AsNumpy(b_cols)

    s_v = np_s[signal_ref_col]
    b_v = np_b[background_ref_col]

    nsig = len(s_v)
    nbkg = len(b_v)

    s_w = np_s[signal_weight_col] if signal_weight_col else np.ones(nsig, dtype=np.float64)
    b_w = np_b[backgound_weight_col] if backgound_weight_col else np.ones(nbkg, dtype=np.float64)

    if nsig == 0 or nbkg == 0:
        raise ValueError(f"While making {title}: \n \tEmpty signal or background arrays")

    if len(s_w) != nsig or len(b_w) != nbkg:
        raise ValueError(f"While making {title}: \n \tMismatching event counts for value and weight arrays")

    descending = -1 if is_higher else 1
    s_v = [np.sort(event)[::descending] for event in s_v]
    b_v = [np.sort(event)[::descending] for event in b_v]

    # ---- Build threshold grid ----
    # Start with all unique combined values; reduce with quantiles if too many.
    thr_all = np.unique(np.concatenate([np.concatenate(s_v), np.concatenate(b_v)]))
    if thr_all.size > ngrid:
        # Quantile-based grid over combined distribution
        qs = np.linspace(0.0, 1.0, ngrid)
        thr_all = np.quantile(thr_all, qs)

    # function to compute signal and background counts given sig, bkg array and threshold
    def compute_counts(sarr, sweight, barr, bweight, thr, is_higher=True):
        op = operator.gt if is_higher else operator.lt
        s_mask = np.array([op(event[0], thr) for event in sarr])
        b_mask = np.array([op(event[0], thr) for event in barr])
        return np.sum(sweight * s_mask), np.sum(bweight * b_mask)
    
    stot = np.sum(s_w)
    btot = np.sum(b_w)

    seff_list = []
    brej_list = []
    for thr in thr_all:
        sc, bc = compute_counts(s_v, s_w, b_v, b_w, thr, is_higher)
        seff_list.append(sc / stot)
        brej_list.append((btot-bc) / btot)
        # print(round(thr, 5), round(sc/stot, 5), round((btot-bc)/btot, 5))

    auc = float(np.trapezoid(seff_list, brej_list))

    # Build a ROOT TGraph (FPR on x-axis, TPR on y-axis)
    seff_arr = np.asarray(seff_list, dtype=np.float64)
    brej_arr = np.asarray(brej_list, dtype=np.float64)
    c = ROOT.TCanvas("c1", "c1", 800, 600)
    g = ROOT.TGraph(seff_arr.size, seff_arr, brej_arr)
    g.SetName(f"{title}")
    g.SetTitle(f"{title}; AUC = {auc}")
    g.SetLineColor(ROOT.kBlue + 2)
    g.SetLineWidth(3)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(ROOT.kBlue + 2)
    g.Draw("ap")
    c.SetLogy()
    c.SaveAs(f"{title}.pdf")

    return

# Example ROC use
# rdf_g.roc_from_rdfs(dfER, dfER, f'{sufElER}_puppiIso', f'{sufElER}_puppiIso', order="lower", ngrid=20)
