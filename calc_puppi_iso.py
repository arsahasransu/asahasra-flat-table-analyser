import numpy as np

from ROOT import RDataFrame


def recalculate_puppi_iso(df: RDataFrame, elcoll: str, puppicoll: str):

    # dRmin_list = np.arange(0.01, 0.2, 0.01)
    dRmin_list = [0.01]
    dRmax = 0.4
    for dRmin in dRmin_list:
        getisostr = f'calcisoannularcone({elcoll}_pt, {elcoll}_eta, {elcoll}_caloEta,\
                        {elcoll}_phi, {elcoll}_caloPhi,\
                        {puppicoll}_pt, {puppicoll}_eta, {puppicoll}_phi, {puppicoll}_pdgId,\
                        {dRmin}, {dRmax})'
        df = df.Define(f'{elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple', getisostr)

        df = df.Define(f'{elcoll}_reisotot_dRmin{str(dRmin).replace('.', '_')}_absPuppiIso',
                       f'std::get<0>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
        df = df.Define(f'{elcoll}_reisochg_dRmin{str(dRmin).replace('.', '_')}_absPuppiIso',
                       f'std::get<7>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
        df = df.Define(f'{elcoll}_reisonut_dRmin{str(dRmin).replace('.', '_')}_absPuppiIso',
                       f'std::get<8>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
        df = df.Define(f'{elcoll}_reisotot_dRmin{str(dRmin).replace('.', '_')}_puppiIso',
                       f'std::get<9>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
        df = df.Define(f'{elcoll}_reisochg_dRmin{str(dRmin).replace('.', '_')}_puppiIso',
                       f'std::get<15>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
        df = df.Define(f'{elcoll}_reisonut_dRmin{str(dRmin).replace('.', '_')}_puppiIso',
                       f'std::get<16>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
    return df


def recalcpuppiiso_comps_oneel(df: RDataFrame, elcoll: str, puppicoll: str):

    # dRmin_list = np.round(np.arange(0.01, 0.21, 0.01), 2)
    dRmin_list = [0.01]
    # dRmax_list = np.round(np.arange(0.1, 0.6, 0.1), 2)
    dRmax_list = [0.4]
    for dRmin in dRmin_list:
        for dRmax in dRmax_list:
            if dRmin >= dRmax:
                continue
            getisostr = f'calcisoanncone_singleobj({elcoll}_pt, {elcoll}_eta, {elcoll}_caloEta,\
                            {elcoll}_phi, {elcoll}_caloPhi,\
                            {puppicoll}_pt, {puppicoll}_eta, {puppicoll}_phi, {puppicoll}_pdgId,\
                            {dRmin}, {dRmax})'
            df = df.Define(f'{elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple', getisostr)

            df = df.Define(f'{elcoll}_reisotot_dRmin{str(dRmin).replace('.', '_')}_absPuppiIso',
                           f'std::get<0>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
            df = df.Define(f'{elcoll}_reisochg_dRmin{str(dRmin).replace('.', '_')}_absPuppiIso',
                           f'std::get<7>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
            df = df.Define(f'{elcoll}_reisonut_dRmin{str(dRmin).replace('.', '_')}_absPuppiIso',
                           f'std::get<8>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
            df = df.Define(f'{elcoll}_reisotot_dRmin{str(dRmin).replace('.', '_')}_puppiIso',
                           f'std::get<9>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
            df = df.Define(f'{elcoll}_reisochg_dRmin{str(dRmin).replace('.', '_')}_puppiIso',
                           f'std::get<15>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')
            df = df.Define(f'{elcoll}_reisonut_dRmin{str(dRmin).replace('.', '_')}_puppiIso',
                           f'std::get<16>({elcoll}_dRmin{str(dRmin).replace('.', '_')}_reisotuple)')

    return df
