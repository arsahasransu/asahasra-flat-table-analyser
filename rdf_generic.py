import os
import pickle
import re


from ROOT import RDataFrame
import ROOT

import numpy as np

from varmetadata import linkvartohist
from pypkg.my_py_generic_utils import time_eval

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


def add_hists_singlecollection(df: RDataFrame, 
                               histograms: list,  # Adds new histograms to this existing list
                               collection: str,   # ID of the collection for variables
                               sreg: str = ''):

    """ Adds histograms from input collection name. Autolists variables in the dataframe.
     If columns don't have the format COLLECTION_VARIABLE, a regular expression can be
     used to specify COLLECTION_REGEX_VARIABLE. Morever the variables should be pre-defined
     in the lnkvartohist list. The function mutates the histograms list.
    """

    # Infer the filter ID of the dataframe
    filtsuf = df.GetFilterNames()[-1]+'_' if df.GetFilterNames() else 'base_'

    # Infer the available variables for the collection and selection ID
    collkey = collection
    dfcolumnnames = list(df.GetColumnNames())

    # Specify reg when a special column is created of the form COLLECTION_RANDOMSTRING_VARNAME
    if sreg != '':
        fullcolumnregex = rf'({collection})_({sreg})_(\w+)'
        selection_regex = re.compile(fullcolumnregex)
        for dfcolname in dfcolumnnames:
            m = re.fullmatch(selection_regex, str(dfcolname))
            if m is None:
                continue
            varname = m.group(3)
            hname = m.group(0)
            if not varname in linkvartohist:
                continue
            hatr = linkvartohist[varname]
            histograms.extend(
                [df.Histo1D((f'{filtsuf}{hname}', hatr[3], hatr[0], hatr[1], hatr[2]), hname)]
            )

    # Default picks columns of the form COLLECTION_VARNAME
    else:
        collkey_vars = [columnname.split('_')[-1] for columnname in dfcolumnnames if columnname.startswith(collkey+'_')]
        collkey_vars = list(set(collkey_vars))
        for var in collkey_vars :
            hname = f'{collkey}_{var}'
            if var not in linkvartohist or hname not in dfcolumnnames:
                continue
            hatr = linkvartohist[var]
            histograms.extend(
                [df.Histo1D((f'{filtsuf}{hname}', hatr[3], hatr[0], hatr[1], hatr[2]), hname)]
            )


def add_hists_multiplecolls(df: RDataFrame, histograms: list, collections: list):
    for collection in collections:
        if ':' in collection:
            collectionparts = collection.split(':')
            add_hists_singlecollection(df, histograms, collectionparts[0], collectionparts[1])
        else:
            add_hists_singlecollection(df, histograms, collection)


@time_eval
def load_rdf_snapshot_from_root(
    rootpath: str,
    treename: str = 'snapshot',
    step_size: int = 100_000
):
    """
    Load ROOT snapshot tree → dict of awkward arrays.

    Works with both uproot 4 (library="ak") and uproot 5 (arrays() already
    returns awkward, no library kwarg — and no interpretation_executor).
    Fixed-entry branches are re-batched per chunk so ak.concatenate can
    rebuild a regular layout across chunks (required for RNTuples, whose
    per-chunk arrays are jagged over cluster boundaries).
    """
    import uproot
    import awkward as ak
    from concurrent.futures import ThreadPoolExecutor

    # uproot 5: RNTuple/LegacyTTree arrays() take neither library nor
    # interpretation_executor; uproot 4: both exist.
    is_uproot5 = int(uproot.__version__.split('.')[0]) >= 5

    executor = ThreadPoolExecutor()   # uproot picks thread count

    try:
        with uproot.open(rootpath) as f:
            tree = f[treename]
            branches = tree.keys()

            kwargs = {
                'decompression_executor': executor,
            }
            if not is_uproot5:
                kwargs['library'] = 'ak'   # use awkward arrays directly
                kwargs['interpretation_executor'] = executor

            # --- accumulate per-branch chunks --------------------------
            accumulators: dict[str, list] = {b: [] for b in branches}

            for batch in tree.iterate(
                branches,
                step_size=step_size,
                **kwargs,
            ):
                for b in branches:
                    arr = batch[b]
                    # Ragged (variable-length) branches concatenate directly.
                    # Fixed-entry branches must be re-batched so each chunk
                    # has the same outer length as its chunk.
                    if hasattr(arr, 'layout') and arr.layout.is_leaf:
                        if ak.num(arr, axis=0) == 1:
                            accumulators[b].append(arr[0])
                        else:
                            accumulators[b].extend(ak.unflatten(arr, ak.num(arr)))
                    else:
                        accumulators[b].append(arr)

            # --- concatenate chunks ------------------------------------
            result = {}
            for b in branches:
                chunks = accumulators[b]
                if not chunks:
                    result[b] = ak.Array([])
                    continue
                result[b] = ak.concatenate(chunks)

    finally:
        executor.shutdown(wait=False)

    return result


def save_rdf_snapshot(df: RDataFrame, cols: list[str], savename: str, *, recreate = False):
    rootfilename    = f'{savename}_snapshot.root'
    treename        = 'snapshot'

    print("Saving: " + ", ".join(cols) + f" to {rootfilename}")

    # --- Stage 1: Snapshot to ROOT (multithreaded C++ side, no Python GIL) ---
    df.Snapshot(treename, rootfilename, cols)
