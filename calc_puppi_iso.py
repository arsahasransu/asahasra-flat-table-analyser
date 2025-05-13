import numpy as np

from ROOT import RDataFrame

def recalculate_puppi_iso(df:RDataFrame, elcoll:str, puppicoll:str):

    dRmin_list = np.arange(0.01, 0.2, 0.01)
    dRmax = 0.4
    for dRmin in dRmin_list:
        getisostr  = f'calcisoannularcone({elcoll}_pt, {elcoll}_eta, {elcoll}_phi,\
                        {puppicoll}_pt, {puppicoll}_eta, {puppicoll}_phi, {dRmin}, {dRmax})'
        df = df.Define(f'{elcoll}_Re_dRmin{str(dRmin).replace('.', '_')}_puppiIso', getisostr)

    return df
