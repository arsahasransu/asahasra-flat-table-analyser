import string
from typing import Dict, List
import warnings

from ROOT import RDataFrame

from varmetadata import linkvartohist


# Filter on a single collections with a list of potential variables
def define_newcollection(df:RDataFrame, collection:string, selection:string, sid:string):
    dfcolumnnames = list(df.GetColumnNames())
    collkey_vars = [columnname.split('_')[-1] for columnname in dfcolumnnames if columnname.startswith(collection)]
    for var in collkey_vars:
        if ( f'{collection}_{var}' in dfcolumnnames ) and ( 'RVec' in df.GetColumnType(f'{collection}_{var}') ):
            df = df.Define(f'{collection}_{sid}_{var}',
                        f'{collection}_{var}[{selection}]')
    
    df = df.Define(f'{collection}_{sid}_n', f'{collection}_{sid}_pt.size()')

    return df


# Add histograms for a collection of variables
# Input: optionally the selection ID of the collection
# Infers the filter ID and the list of available variables for the collection
def add_hists_singlecollection(df:RDataFrame, histograms:List, collection:string, sid:string=''):

    # Infer the filter ID of the dataframe
    fid = ''
    # TODO LINE BELOW GIVES MORE TIME DELAY THAN EXPECTED - REASON UNKNOWN
    filternameslist = df.GetFilterNames()
    if len(filternameslist) != 0:
        fid = filternameslist[-1] + '_'

    # Infer the available variables for the collection and selection ID
    collkey = collection if sid=='' else f'{collection}_{sid}'
    dfcolumnnames = list(df.GetColumnNames())
    collkey_vars = [columnname.split('_')[-1] for columnname in dfcolumnnames if columnname.startswith(collkey)]
    histograms.extend([df.Histo1D((f'{fid}{collkey}_{var}', linkvartohist[var][3],
                                   linkvartohist[var][0], linkvartohist[var][1], linkvartohist[var][2]),
                                  f'{collkey}_{var}') for var in collkey_vars if var in linkvartohist])
    

def add_hists_multiplecolls(df:RDataFrame, histograms:List, collections:List):
    for collection in collections:
        if ':' in collection:
            collectionparts = collection.split(':')
            add_hists_singlecollection(df, histograms, collectionparts[0], collectionparts[1])
        else:
            add_hists_singlecollection(df, histograms, collection)
