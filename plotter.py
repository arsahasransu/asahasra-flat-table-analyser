import glob
import os
import pathlib
import re
import yaml

import ROOT
from ROOT import TH1, TH1D
from ROOT import TCanvas, TEfficiency, TFile

from an_specific_utilities import make_plotnorm_by_scheme
from my_py_generic_utils import recreate_dir
from plotBeautifier import makePngPlot

# Enable batch mode
ROOT.gROOT.SetBatch(True)


def make_effs(dirname:str, base:TH1, target:TH1, var:str='pt'):
    """
    Make efficiency plots from base and target TH1 ROOT histograms
    """

    eff = TEfficiency(target, base)

    c = TCanvas("c", "c", 800, 600)
    eff.Draw()

    recreate_dir(f"plots/{dirname}/eff")
    c.SaveAs(f"plots/{dirname}/eff/{var}_eff.png")


def plotter():
    """
    The class to handle plotting for all different settings and possibilities
    in the analysis.
    * Compare same selection variables between different samples - DEFAULT
    * Compare variables for the same sample with different selection
    * Currently plots TH1 variables only
    * Possible normalisation schemes - default_no_norm, hist_integral, summed_components
    """

    with open('./plot_config.yaml') as config_file:
        config = yaml.safe_load(config_file)

    # Load the general settings
    fdir = config['general']['hist_folder']  # Path to histogram directory
    samples = config['general']['samples']  # Dictionary of sample names and file paths

    for plotdir, plotconfig in config.items():
        if plotdir == 'general':
            continue

        try:
            recreate_dir(f"plots/{plotdir}")
        except Exception as e:
            print(f"Exception {e} occurred while creating plot directory {plotdir}")

        if "norm_scheme" in plotconfig.keys():
            norm_scheme = plotconfig['norm_scheme']
        else:
            norm_scheme = {'method': 'default_no_norm'}

        if 'type' in list(plotconfig.keys()) and plotconfig['type'] == 'same_sample_different_collection':
            # Use plotconfig['hist_tag'] to open the right collection name and gather attributes to plot
            pfxs = plotconfig['hist_tag']['prefix']
            colls = plotconfig['hist_tag']['collection']
            if norm_scheme['method'] == "summed_components":
                raise NotImplementedError("Unsupported normalisation scheme \"summed_components\" for "\
                                          f"plotting type \"same_sample_different_collection\" in {plotdir}")

            for folder in plotconfig['ftags']:
                fname = samples[folder]
                recreate_dir(f"plots/{plotdir}/{folder}")
                file = ROOT.TFile.Open(f'{fdir}/{fname}')
                for coll in colls:
                    plotleg = []
                    collnames = []
                    for pfxk, pfxv in pfxs.items():
                        plotleg.append(pfxk)
                        collnames.append(coll.format(tag=pfxv))
                    varnames = list(set([hkey.GetName().split('_')[-1] for hkey in file.GetListOfKeys()if hkey.GetName().startswith(collnames[0]+'_')]))
                    allhists = []
                    for varname in varnames:
                        hists = []
                        for collname in collnames:
                            hist = file.Get(f'{collname}_{varname}')
                            # Check if TH1 is in name of the histogram
                            if hist and hist.InheritsFrom(ROOT.TH1.Class()):
                                hists.append(hist)
                        if len(plotleg) == len(hists):
                            allhists.append(hists)

                    normcnts = make_plotnorm_by_scheme(allhists, norm_scheme['method'])
                    
                    for hist, normcnt in zip(allhists, normcnts):
                        makePngPlot(hist, f'plots/{plotdir}/{folder}', 'autoCompPlot', plotleg, normcnt)
                file.Close()

        # DEFAULT - Plot same/analogous collection for different samples            
        else:
            files = []  # Each row is a different sample
            legend = []  # Same number of entries as samples
            fcollnames = []  # Each collumn is a different collection
            sel_col_re = re.compile(r"^([A-Za-z0-9_]+)\[([A-Za-z0-9_, ]+)\]$")
            alphanum_re = re.compile(r"^[A-Za-z0-9_]+$")
            summed_component = norm_scheme['norm_config']['refsample'] if norm_scheme['method'] == "summed_components" else -1
            summed_component_pos = -1

            for i, (ftag, colls) in enumerate(plotconfig['hist_tag'].items()):
                if summed_component == ftag:
                    summed_component_pos = i
                legend.append(ftag)
                files.append(ROOT.TFile.Open(f'{fdir}/{samples[ftag]}'))

                # Convert pattern to individual colls
                # a_[b, c] => a_b, a_c
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
            if norm_scheme['method'] == "summed_components" and summed_component_pos < 0:
                raise RuntimeError('Summed component index not found in sample list to perform correct normalisation')

            for collnames in fcollnames:
                allkeys = [hkey.GetName() for hkey in files[0].GetListOfKeys()]
                vars = list(set([hkey.split('_')[-1] for hkey in allkeys if hkey.startswith(collnames[0]+'_')]))

                groupbyvars_hists = []
                for var in vars:
                    groupbyvars_hist = []
                    for file, collname in zip(files, collnames):
                        hist = file.Get(f'{collname}_{var}')
                        # Check if TH1 is in name of the histogram
                        if hist and hist.InheritsFrom(ROOT.TH1.Class()):
                            groupbyvars_hist.append(hist)
                    if len(legend) == len(groupbyvars_hist):
                        groupbyvars_hists.append(groupbyvars_hist)

                if(norm_scheme['method'] == "summed_components"):
                    normcnts = make_plotnorm_by_scheme(groupbyvars_hists, norm_scheme['method'],
                                                       summed_sample_pos=summed_component_pos, byevent_vartag='n')
                else:
                    normcnts = make_plotnorm_by_scheme(groupbyvars_hists, norm_scheme['method'])

                for groupbyvar_hist,normcnt in zip(groupbyvars_hists,normcnts):
                    makePngPlot(groupbyvar_hist, f'plots/{plotdir}', 'autoCompPlot', legend, normcnt)


def auto_singlehist_plotter():

    histdirs = glob.glob("./OutHistoFiles/hists_*.root")

    for file_n in histdirs:
        file = TFile.Open(file_n)
        basefilename = os.path.basename(file_n)
        file_prefix = basefilename[:-5]

        # Remake the histogram autoplot directory
        try:
            recreate_dir(f'plots/{file_prefix}/autoplots')
        except Exception as e:
            print(f"Exception {e} occurred while creating plot directory {file_prefix}/autoplots")

        for histKey in file.GetListOfKeys():
            histObj = histKey.ReadObj()
            if isinstance(histObj, ROOT.TH1D):
                histTh1Obj = TH1D(histKey.ReadObj())

                makePngPlot([histTh1Obj], f'plots/{file_prefix}/autoplots', 'autoSinglePlot')


if __name__ == "__main__":
    plotter()
    # auto_singlehist_plotter()
