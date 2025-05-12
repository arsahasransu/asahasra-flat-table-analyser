import os
import glob

import ROOT
from ROOT import TCanvas, TEfficiency, TH1D, TFile

from plotBeautifier import makePngPlot

# Enable batch mode
ROOT.gROOT.SetBatch(True)

def make_effs(dirname, base, target, var='pt'):

    eff = TEfficiency(target, base)

    c = TCanvas("c", "c", 800, 600)
    eff.Draw()
    if not os.path.exists(f"{dirname}/eff"):
        os.makedirs(f"{dirname}/eff")
    c.SaveAs(f"{dirname}/eff/{var}_eff.png")


def plotter():

    file_list = os.listdir("./OutHistoFiles")
    
    for fname in file_list:
        file = TFile.Open(f"./OutHistoFiles/{fname}")

        base_h = file.Get('GenElPromptEB_pt')
        target_h = file.Get('GenElPromptEB_PuppiElEBmatched_pt')
        dirname = fname.split('.')[0]
        make_effs(dirname, base_h, target_h)

        base_h = file.Get('GenElPromptEB_eta')
        target_h = file.Get('GenElPromptEB_PuppiElEBmatched_eta')
        dirname = fname.split('.')[0]
        make_effs(dirname, base_h, target_h, 'eta')

        base_h = file.Get('GenElPromptEB_phi')
        target_h = file.Get('GenElPromptEB_PuppiElEBmatched_phi')
        dirname = fname.split('.')[0]
        make_effs(dirname, base_h, target_h, 'phi')


def auto_plotter():

    histdirs = glob.glob("./OutHistoFiles/hists_*.root")


    for file_n in histdirs:
        file = TFile.Open(file_n)
        basefilename = os.path.basename(file_n)
        file_prefix = basefilename[:-5]

        for histKey in file.GetListOfKeys():
            histObj = histKey.ReadObj()
            if isinstance(histObj, ROOT.TH1D):
                histTh1Obj = TH1D(histKey.ReadObj())

                makePngPlot(histTh1Obj, f'{file_prefix}/autoplots', 'autoSinglePlot')


if __name__ == "__main__":
    # plotter()
    auto_plotter()