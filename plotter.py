import copy
import glob
import os
import pathlib
import shutil
import yaml

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

    with open('./plot_config.yaml') as config_file:
        config = yaml.safe_load(config_file)

    fdir = config['general']['hist_folder']
    for k, v in config.items():
        if k == 'general':
            continue

        # Use k to create a directory
        pdir = pathlib.Path(f'{k}')
        if pdir.exists() and pdir.is_dir():
            shutil.rmtree(pdir)
        pdir.mkdir(parents=True, exist_ok=True)

        # Use v['samples'] to open the right root file
        legend = list(v['samples'].keys())
        fnames = v['samples'].values()
        files = [ROOT.TFile.Open(f'{fdir}/{fname}') for fname in fnames]

        # Use v['hist_tag'] to open the right collection name and gather attributes to plot
        pfxs = v['hist_tag']['prefix']
        colls = v['hist_tag']['collection']

        if v['type'] == 'same_sample_different_collection':
            for fname, folder in zip(fnames, legend):
                pathlib.Path(k+'/'+folder).mkdir(parents=True, exist_ok=True)
                file = ROOT.TFile.Open(f'{fdir}/{fname}')
                for coll in colls:
                    plotleg = []
                    collnames = []
                    for pfxk, pfxv in pfxs.items():
                        plotleg.append(pfxk)
                        collnames.append(coll.format(tag=pfxv))
                    varnames = [hkey.GetName().split('_')[-1] for hkey in file.GetListOfKeys()if hkey.GetName().startswith(collnames[0]+'_')]
                    for varname in varnames:
                        hnames = [collname+'_'+varname for collname in collnames]
                        hobjs = [file.Get(hname) for hname in hnames]

                        makePngPlot(hobjs, f'{k}/{folder}', 'autoCompPlot', plotleg)
                file.Close()

        # The default behaviour is to plot for same collection different samples            
        else:
            for coll in colls:
                file0_hname = f'{pfxs[legend[0]]}_{coll}'
                exactcolls = ['_'.join(hkey.GetName().split('_')[:-1]) for hkey in files[0].GetListOfKeys() if hkey.GetName().startswith(file0_hname)]
                exactcolls = list(set(exactcolls))
                for exactcoll in exactcolls:
                    varnames = [hkey.GetName().split('_')[-1] for hkey in files[0].GetListOfKeys() if hkey.GetName().startswith(exactcoll+'_')]
                    for varname in varnames:
                        collname = '_'.join(exactcoll.split('_')[1:])
                        hnames = [f'{pfxs[leg]}_{collname}_{varname}' for leg in legend]
                        hobjs = [file.Get(hname) for file, hname in zip(files, hnames)]

                        makePngPlot(hobjs, f'{k}', 'autoCompPlot', legend)

        # for row in v['hist_tag'][1:]:
        #     colls_copied = copy.deepcopy(colls)
        #     colls.clear()
        #     for colsub in row:
        #         for col in colls_copied:
        #             colls.append(f'{col}_{colsub}')

        # for coll in colls:
        #     hnames = [hkey.GetName() for hkey in files[0].GetListOfKeys() if hkey.GetName().startswith(coll)]
        #     for hname in hnames:
        #         hobjs = [file.Get(hname) for file in files]

        #         makePngPlot(hobjs, f'{k}', 'autoCompPlot', legend)

    # Legacy Efficiency Plot Maker

    # for fname in file_list:
    #     file = TFile.Open(f"./OutHistoFiles/{fname}")

    #     base_h = file.Get('GenElPromptEB_pt')
    #     target_h = file.Get('GenElPromptEB_PuppiElEBmatched_pt')
    #     dirname = fname.split('.')[0]
    #     make_effs(dirname, base_h, target_h)

    #     base_h = file.Get('GenElPromptEB_eta')
    #     target_h = file.Get('GenElPromptEB_PuppiElEBmatched_eta')
    #     dirname = fname.split('.')[0]
    #     make_effs(dirname, base_h, target_h, 'eta')

    #     base_h = file.Get('GenElPromptEB_phi')
    #     target_h = file.Get('GenElPromptEB_PuppiElEBmatched_phi')
    #     dirname = fname.split('.')[0]
    #     make_effs(dirname, base_h, target_h, 'phi')


def auto_singlehist_plotter():

    histdirs = glob.glob("./OutHistoFiles/hists_*.root")

    for file_n in histdirs:
        file = TFile.Open(file_n)
        basefilename = os.path.basename(file_n)
        file_prefix = basefilename[:-5]

        # Remake the histogram autoplot directory
        outdir = pathlib.Path(f'{file_prefix}/autoplots')
        if outdir.exists() and outdir.is_dir():
            shutil.rmtree(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        for histKey in file.GetListOfKeys():
            histObj = histKey.ReadObj()
            if isinstance(histObj, ROOT.TH1D):
                histTh1Obj = TH1D(histKey.ReadObj())

                makePngPlot([histTh1Obj], f'{file_prefix}/autoplots', 'autoSinglePlot')


if __name__ == "__main__":
    plotter()
    # auto_singlehist_plotter()
