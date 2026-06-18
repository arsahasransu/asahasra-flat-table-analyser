import os
import pickle
import re


from ROOT import RDataFrame
import ROOT

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from concurrent.futures import ThreadPoolExecutor

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


def save_rdf_snapshot_to_parquet(df: RDataFrame, cols: list[str], savename: str, *, recreate = False):
    """Save RDataFrame columns to .parquet — RVec<float> → pyarrow list<item: float32>
    
    Columns that are RVec (object dtype) → variable-length list per event.
    Scalar columns preserved as-is.
    """
    parquetfilename = f'{savename}.parquet'
    print("Saving: " + ", ".join(cols) + f" to {parquetfilename}")

    np_dict = df.AsNumpy(cols)
    table_dict = {}

    def convert_column(item):
        k, arr = item
        if arr.dtype == object:
            return k, pa.array(arr.tolist(), type=pa.list_(pa.float32()))
        return k, pa.array(arr)

    with ThreadPoolExecutor() as executor:
        results = executor.map(convert_column, np_dict.items())

    table_dict = dict(results)
    # for k, arr in np_dict.items():
    #     if arr.dtype == object:
    #         table_dict[k] = pa.array(arr.tolist(), type=pa.list_(pa.float32()))
    #     else:
    #         table_dict[k] = pa.array(arr)

    table = pa.table(table_dict)

    if os.path.exists(parquetfilename) and not recreate:
        existing = pq.read_table(parquetfilename)
        schema_names = set(existing.schema.names)

        new_cols = {}
        for k, arr in table_dict.items():
            if k in schema_names:
                raise RuntimeError(f'Key {k} exists in {parquetfilename}')
            new_cols[k] = arr

        # Single table construction, not iterative append
        merged = existing.append_column  # wrong pattern — do this instead:
        merged = pa.table({**{n: existing.column(n) for n in existing.schema.names}, **new_cols})
        pq.write_table(merged, parquetfilename)        
    #     existing = pq.read_table(parquetfilename)
    #     existing_names = set(existing.column_names)
    #     for k, arr in table_dict.items():
    #         if k not in existing_names:
    #             existing = existing.append_column(pa.field(k, type=arr.type), arr)
    #         else:
    #             raise RuntimeError(f'Key {k} exists in file {parquetfilename}')
    #     pq.write_table(existing, parquetfilename)
    # else:
    #     pq.write_table(table, parquetfilename)


def load_rdf_snapshot_from_parquet(parquetpath: str) -> dict[str, np.ndarray]:
    """Load .parquet snapshot → dict of numpy arrays.
    list<float> columns → object arrays with per-event vectors preserved."""
    table = pq.read_table(parquetpath)
    result = {}
    for col in table.column_names:
        arr = table.column(col)
        if pa.types.is_list(arr.type):
            result[col] = np.array(arr.to_pylist(), dtype=object)
        else:
            result[col] = arr.to_numpy()
    return result


def save_rdf_snapshot(df: RDataFrame, cols: list[str], savename: str, *, recreate = False):
    parquetfilename = f'{savename}.parquet'
    rootfilename    = f'{savename}_snapshot.root'
    treename        = 'snapshot'

    print("Saving: " + ", ".join(cols) + f" to {parquetfilename}")

    # --- Stage 1: Snapshot to ROOT (multithreaded C++ side, no Python GIL) ---
    df.Snapshot(treename, rootfilename, cols)
