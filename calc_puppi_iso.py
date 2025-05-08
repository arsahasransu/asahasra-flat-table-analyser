import numpy as np


def get_iso(df, elcoll, puppicoll):

    dRmin_list = np.arange(0.01, 0.2, 0.01)
    dRmax = 0.4
    for dRmin in dRmin_list:
        getisostr  = f'calcpuppiiso({elcoll}_pt, {elcoll}_eta, {elcoll}_phi,\
                        {puppicoll}_pt, {puppicoll}_eta, {puppicoll}_phi, {dRmin}, {dRmax})'
        df = df.Define(f'{elcoll}_recalcPuppiIso_dRmin{str(dRmin).replace('.', '_')}', getisostr)

    return df