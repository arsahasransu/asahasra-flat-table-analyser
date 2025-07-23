import numpy as np

from ROOT import RDataFrame


def recalculate_puppi_iso(df: RDataFrame, elcoll: str, puppicoll: str):

    dRmin_list = np.arange(0.01, 0.2, 0.01)
    dRmax = 0.4
    for dRmin in dRmin_list:
        getisostr = f'calcisoannularcone({elcoll}_pt, {elcoll}_eta, {elcoll}_phi,\
                        {puppicoll}_pt, {puppicoll}_eta, {puppicoll}_phi, {dRmin}, {dRmax})'
        df = df.Define(f'{elcoll}_Re_dRmin{str(dRmin).replace('.', '_')}_puppiIso', getisostr)

    return df


def recalcpuppiiso_comps_oneel(df: RDataFrame, elcoll: str, puppicoll: str):

    dRmin_list = np.round(np.arange(0.01, 0.21, 0.01), 2)
    dRmax_list = np.round(np.arange(0.1, 0.6, 0.1), 2)
    for dRmin in dRmin_list:
        for dRmax in dRmax_list:
            if dRmin >= dRmax:
                continue
            getisostr = f'calcisoanncone_singleobj({elcoll}_pt, {elcoll}_eta, {elcoll}_caloEta,\
                            {elcoll}_phi, {elcoll}_caloPhi,\
                            {puppicoll}_pt, {puppicoll}_eta, {puppicoll}_phi, {puppicoll}_pdgId,\
                            {dRmin}, {dRmax})'
            for i, isostr in enumerate(['isotot', 'iso11', 'iso13', 'iso22', 'iso130', 'iso211', 'isooth']):
                collname = f'{elcoll}_Re_'
                collname += f'dRmin{str(dRmin).replace('.', '_')}_'
                collname += f'dRmax{str(dRmax).replace('.', '_')}_'
                df = df.Define(f'{collname}{isostr}_puppiIso', f'std::get<{i}>({getisostr})')

    return df
