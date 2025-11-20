import glob
import os
import pathlib
import re
import yaml

import ROOT
from ROOT import TCanvas, TEfficiency, TH1D, TFile

from an_specific_utilities import make_plotnorm_by_scheme
from my_py_generic_utils import recreate_dir
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
    samples = config['general']['samples']
    norm_scheme = {'method': 'default_no_norm'}

    for plotdir, plotconfig in config.items():
        if plotdir == 'general':
            continue

        try:
            recreate_dir(plotdir)
        except Exception as e:
            print(f"Exception {e} occurred while creating plot directory {plotdir}")

        if "norm_scheme" in plotconfig.keys():
            norm_scheme = plotconfig['norm_scheme']
        else:
            norm_scheme = {'method': 'default_no_norm'}

        # Use plotconfig['hist_tag'] to open the right collection name and gather attributes to plot
        if 'type' in list(plotconfig.keys()) and plotconfig['type'] == 'same_sample_different_collection':
            pfxs = plotconfig['hist_tag']['prefix']
            colls = plotconfig['hist_tag']['collection']

            for folder in plotconfig['ftags']:
                fname = samples[folder]
                pathlib.Path(plotdir+'/'+folder).mkdir(parents=True, exist_ok=True)
                file = ROOT.TFile.Open(f'{fdir}/{fname}')
                for coll in colls:
                    plotleg = []
                    collnames = []
                    for pfxk, pfxv in pfxs.items():
                        plotleg.append(pfxk)
                        collnames.append(coll.format(tag=pfxv))
                    varnames = list(set([hkey.GetName().split('_')[-1] for hkey in file.GetListOfKeys()if hkey.GetName().startswith(collnames[0]+'_')]))
                    for varname in varnames:
                        hnames = [collname+'_'+varname for collname in collnames]
                        hobjs = [file.Get(hname) for hname in hnames]

                        normcnt = make_plotnorm_by_scheme(hobjs, norm_scheme['method'])
                        makePngPlot(hobjs, f'{plotdir}/{folder}', 'autoCompPlot', plotleg, normcnt)
                file.Close()

        # The default behaviour is to plot for same/analogous collection different samples            
        else:
            files = []
            fcollnames = []
            legend = []
            sel_col_re = re.compile(r"^([A-Za-z0-9_]+)\[([A-Za-z0-9_, ]+)\]$")
            alphanum_re = re.compile(r"^[A-Za-z0-9_]+$")

            for ftag, colls in plotconfig['hist_tag'].items():
                legend.append(ftag)
                files.append(ROOT.TFile.Open(f'{fdir}/{samples[ftag]}'))
                # Make the individual collection names from the composite input
                getselmatch = sel_col_re.match(colls)
                if getselmatch:
                    sel = getselmatch.group(1)
                    allcolls = map(str.strip, getselmatch.group(2).split(','))
                    for i, coll in enumerate(allcolls):
                        if len(fcollnames) <= i:
                            fcollnames.append([])
                        fcollnames[i].append(f"{sel}{coll}")
                elif alphanum_re.match(colls):
                    if len(fcollnames) == 0:
                        fcollnames.append([])
                    fcollnames[0].append(colls)
                else:
                    raise RuntimeError("Bad collection name format: ", colls)
            fcollnames_rowlength = [len(row) for row in fcollnames]
            if len(set(fcollnames_rowlength)) > 1:
                raise RuntimeError("Inconsistent number of collections given to plot between the ftags: ",\
                                   fcollnames_rowlength)
            for collnames in fcollnames:
                allkeys = [hkey.GetName() for hkey in files[0].GetListOfKeys()]
                vars = list(set([hkey.split('_')[-1] for hkey in allkeys if hkey.startswith(collnames[0]+'_')]))

                groupbyvars_hists = [[file.Get(f'{collname}_{var}') for file,collname in zip(files,collnames)] for var in vars]
                print(len(groupbyvars_hists), len(groupbyvars_hists[0]))
                if(norm_scheme['method'] == "summed_components"):
                    normcnts = make_plotnorm_by_scheme(groupbyvars_hists, norm_scheme['method'],
                                                       summed_sample_pos=2, byevent_vartag='n')
                else:
                    normcnts = make_plotnorm_by_scheme(groupbyvars_hists, norm_scheme['method'])

                for groupbyvar_hist,normcnt in zip(groupbyvars_hists,normcnts):
                    print(normcnt)
                    makePngPlot(groupbyvar_hist, f'{plotdir}', 'autoCompPlot', legend, normcnt)

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
        try:
            recreate_dir(f'{file_prefix}/autoplots')
        except Exception as e:
            print(f"Exception {e} occurred while creating plot directory {file_prefix}/autoplots")

        for histKey in file.GetListOfKeys():
            histObj = histKey.ReadObj()
            if isinstance(histObj, ROOT.TH1D):
                histTh1Obj = TH1D(histKey.ReadObj())

                makePngPlot([histTh1Obj], f'{file_prefix}/autoplots', 'autoSinglePlot')


if __name__ == "__main__":
    plotter()
    # auto_singlehist_plotter()
