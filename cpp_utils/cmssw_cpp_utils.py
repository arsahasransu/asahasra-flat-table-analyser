"""
The utilities as used in CMSSW
"""

import ROOT


def cmssw_cpp_utils():

    STRCPPFUNC_CMSSW_CONSTS = """
        constexpr int INTPHI_PI = 720;
        constexpr float ETAPHI_LSB = M_PI / INTPHI_PI;
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_CMSSW_CONSTS)


    STRCPP_CMSSW_INLINE_FUNCS = """
        inline int makeDR2FromFloatDR(float dr) { return ceil(dr * dr / ETAPHI_LSB / ETAPHI_LSB); }
    """

    ROOT.gInterpreter.Declare(STRCPP_CMSSW_INLINE_FUNCS)
