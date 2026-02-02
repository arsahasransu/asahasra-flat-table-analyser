import os
import pickle
import re


from ROOT import RDataFrame
import ROOT

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


def save_rdf_snapshot_to_pkl(df: RDataFrame, cols: list[str], savename: str, *, recreate = False):
    pklfilename = f'{savename}.pkl'
    print("Saving: " + ", ".join(cols) + f" to {pklfilename}")

    np_dict = df.AsNumpy(cols)
    if os.path.exists(pklfilename) and not recreate:
        with open(pklfilename, "rb") as repklf:
            pkld = pickle.load(repklf)
        for k, v in np_dict.items():
            if k not in pkld:
                pkld[k] = v
            else:
                raise RuntimeError(f'Key {k} exists in file {pklfilename}')
            
        with open(pklfilename, "wb") as pklf:
            pickle.dump(pkld, pklf, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(pklfilename, "wb") as pklf:
            pickle.dump(np_dict, pklf, protocol=pickle.HIGHEST_PROTOCOL)
