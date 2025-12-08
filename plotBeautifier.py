import math as m
import numpy as np
import os
import re
import warnings

import ROOT
from ROOT import TH1

from an_specific_utilities import customisation_conds as cc
from an_specific_utilities import conditionally_modify_plots


# TODO INCOMPLETE
# def adjustBinContentAuto(histo: TH1):
#     # For bins with low stats - adjust bin content to have atleast
#     # minstat percentage of events in the next bin comapred to the previous bin
#     #
#     # Skip bins with zero entries

#     minstat = 0.01*10

#     # Don't change the underflow or overflowbin values
#     binlowedge = [histo.GetBinLowEdge(bin) for bin in range(1, histo.GetNbinsX()+2)]
#     bincontent = [histo.GetBinContent(bin) for bin in range(1, histo.GetNbinsX()+2)]
#     lastgoodbin = 0
#     bincsum = 0.0
#     for nbin, binc in enumerate(bincontent):
#         if nbin == 0 or binc == 0:
#             lastgoodbin = nbin
#             bincsum = 0.0
#             continue
#         if binc > minstat*bincontent[nbin-1]:
#             lastgoodbin = nbin
#             bincsum = 0.0
#             continue
#         newbinc = bincsum + binc
#         if newbinc > minstat*bincontent[nbin-1]:
#     print(len(binlowedge))

#     return histo


def makePlot_isolation(histo: TH1):
    # Generate 100 xbins upto precision 2 in log-scale
    xlow = 0.001
    xhigh = 100
    nbins = 200
    step = (m.log10(xhigh) - m.log10(xlow))/nbins
    bins = [np.round(10**i, 3) for i in np.arange(m.log10(xlow), m.log10(xhigh), step)]
    bins = list(np.sort(list(set(bins))))
    bins.append(xhigh)
    bins_array = np.array(bins, dtype='float64')
    old_bw = histo.GetBinWidth(1)

    # Set isolation value 0 to the first bin to enable log x-axis
    histo.SetBinContent(2, histo.GetBinContent(0) + histo.GetBinContent(1) + histo.GetBinContent(2))
    histo.SetBinContent(1, 0)
    histo.SetBinContent(0, 0)

    histo_rebinned = histo.Rebin(len(bins)-1, histo.GetName()+'_rebinned', bins_array)
    # Combine overflow into last bin
    lbin = histo_rebinned.GetNbinsX()
    histo_rebinned.SetBinContent(lbin, histo_rebinned.GetBinContent(lbin)+histo_rebinned.GetBinContent(lbin+1))
    histo_rebinned.SetBinError(lbin, histo_rebinned.GetBinError(lbin)+histo_rebinned.GetBinError(lbin+1))
    for bin in range(0, histo_rebinned.GetNbinsX()+1):
        binc = histo_rebinned.GetBinContent(bin)
        bine = histo_rebinned.GetBinError(bin)
        bw = histo_rebinned.GetBinWidth(bin)
        ncombined_bins = bw/old_bw
        histo_rebinned.SetBinContent(bin, binc/ncombined_bins)
        histo_rebinned.SetBinError(bin, bine/ncombined_bins)

    histo_rebinned.GetYaxis().SetTitle('Events / bin width')
    return histo_rebinned


def defineXrangeAuto(histo: TH1, includeUnderflow: bool = True, includeOverflow: bool = True):
    # Decide if underflow and overflow bins should be included in the plot
    includeUnderflow = histo.GetBinContent(0) > 0.01*histo.GetBinContent(1)
    includeOverflow = histo.GetBinContent(histo.GetNbinsX()+1) > 0.01*histo.GetBinContent(histo.GetNbinsX()+1)

    # Loop over the bin content and find the fist bin with data in it
    beginBin = 0 if includeUnderflow else 1
    for i in range(beginBin, histo.GetNbinsX()+2):
        if histo.GetBinContent(i) != 0:
            beginBin = i - 1 if i > 0 else i
            break

    # Loop over the bin content and find the last bin with data in it
    endBin = histo.GetNbinsX()+1 if includeOverflow else histo.GetNbinsX()
    for i in range(endBin, beginBin-1, -1):
        if histo.GetBinContent(i) != 0:
            endBin = i + 1 if i < histo.GetNbinsX()+1 else i
            break

    histo.GetXaxis().SetRange(beginBin, endBin)

    return {'uFlow': includeUnderflow, 'oFlow': includeOverflow}


def makeBinLabel(histo: TH1, binNumber: int, label: str):
    x = histo.GetXaxis().GetBinLowEdge(binNumber)
    y = histo.GetBinContent(binNumber)
    labelobj = ROOT.TText(x, y, label)
    return labelobj


def makeGoodSinglePlotAuto(histo: TH1):
    if 'iso' in histo.GetTitle():
        histo = makePlot_isolation(histo)

    objectList = [histo]

    # The optimal x-range is defined as the first and the last bin
    # with non-zero entries including underflow and overflow
    binDecisions = defineXrangeAuto(histo)

    histo.SetLineWidth(4)
    if (binDecisions['uFlow']):
        ulabel = makeBinLabel(histo, 0, 'UFLOW')
        ulabel.SetTextAngle(90)
        ulabel.SetTextAlign(13)
        objectList.append(ulabel)
    if (binDecisions['oFlow']):
        olabel = makeBinLabel(histo, histo.GetNbinsX()+1, 'OFLOW')
        olabel.SetTextAngle(90)
        olabel.SetTextAlign(13)
        objectList.append(olabel)

    histo.GetXaxis().SetTitle(histo.GetTitle())
    histo.SetTitle('')

    return objectList


def histname_to_folder(name):
    # Get the folder and histogram details from the histogram name
    # Assuming name has been saved as DATAFRAME_COLLECTION_<OBJSELECTION_>VAR
    # The OBJSELECTION could be optional
    histMetaNameSplit = re.fullmatch('([A-Za-z0-9]+)_([A-Za-z0-9]+)_([A-Za-z0-9_]+)', name)
    histrdf = histMetaNameSplit.group(1)
    histcoll = histMetaNameSplit.group(2)
    histname = histMetaNameSplit.group(3)

    return histrdf, histcoll, histname


def auto_ylog_decision(hist):
    logy = False

    # If the minimum non-zero bin entry is less than tolerance percentage
    # of the maximum bin entry, plot in log scale
    tolerance = 15
    binentries = [hist.GetBinContent(i) for i in range(hist.GetNbinsX())]
    nonzerobinentries = [binentry for binentry in binentries if binentry != 0]
    binentrymax = max(nonzerobinentries) if len(nonzerobinentries)!=0 else hist.GetNbinsX()+1
    binentrymin = min(nonzerobinentries) if len(nonzerobinentries)!=0 else 0
    if binentrymin < tolerance * 0.01 * binentrymax:
        logy = True

    if len(nonzerobinentries) < 10:
        drawopt = 'TEXT HIST'
    else:
        drawopt = ''

    return logy, drawopt


def autoModifyHists(histlist, histref):
    del histlist[0]
    newhistlist = []
    for hist in histlist:
        if 'iso' in hist.GetTitle():
            newhist = makePlot_isolation(hist)
            newhistlist.append(newhist)
            continue

        nbins = histref.GetNbinsX()
        axis = histref.GetXaxis()
        edges = axis.GetXbins()
        if edges.GetSize() > 0:
            # variable binning
            bin_edges = np.array([edges.At(i) for i in range(edges.GetSize())], dtype='float64')
        else:
            # uniform binning
            xmin = axis.GetXmin()
            xmax = axis.GetXmax()
            bin_edges = np.array([xmin + i * (xmax - xmin) / nbins for i in range(nbins + 1)], dtype='float64')

        hist_rebinned = hist.Rebin(len(bin_edges)-1, hist.GetName(), bin_edges)
        hist_rebinned.GetXaxis().SetRange(histref.GetXaxis().GetFirst(), histref.GetXaxis().GetLast())
        hist_rebinned.SetLineWidth(4)
        hist_rebinned.GetXaxis().SetTitle(hist_rebinned.GetTitle())
        hist_rebinned.SetTitle('')

        newhistlist.append(hist_rebinned)

    return newhistlist


def generateColorPalette(ncolours):

    ROOT.gStyle.SetPalette(ROOT.kCubehelix)

    colours = [ROOT.gStyle.GetColorPalette(int((i+1) * 255/ncolours)) for i in range(ncolours)]
    return colours


def makePngPlot(histList, outputDir: str, plotkey: str, legList=[], normlist=[]):

    histList = conditionally_modify_plots(histList)
    if len(normlist) == 0:
        normlist = [1.0]*len(histList)
    for normv,histv in zip(normlist, histList):
        if normv == 0.0:
            warnings.warn(f"Found null normalisation. Skipping drawing histogram: {histv.GetName()}")

    try:
        if histList[0].Integral() == 0:
            warnings.warn(
                f"Null integral for (atleast) first histogram in list. Skipping drawing histogram: {histList[0].GetName()}")
            return
    except AttributeError as e:
        raise f"{histList[0].GetName()} is of type {histList[0].ClassName()} not TH1"

    if plotkey == 'autoSinglePlot':
        histObj = histList[0]
        histMetaName = histObj.GetName()
        histCopy = histObj.Clone()

        histrdf, histcoll, histname = histname_to_folder(histMetaName)
        drawableObjectList = makeGoodSinglePlotAuto(histCopy)

        c1 = ROOT.TCanvas('c1', 'c1', 800, 600)
        ROOT.gStyle.SetOptStat(0)

        logy, drawopt = auto_ylog_decision(drawableObjectList[0])
        c1.SetLogy(logy)
        if logy:
            drawableObjectList[0].SetMinimum(0.9)
        drawableObjectList[0].Draw(drawopt)

        for i in range(1, len(drawableObjectList)):
            drawableObjectList[i].Draw()

        for canvconds in cc['canvas']:
            res = canvconds[0](histname)
            if res:
                canvconds[1](c1)

        if not os.path.exists(f'{outputDir}/{histrdf}/{histcoll}'):
            os.makedirs(f'{outputDir}/{histrdf}/{histcoll}')
        c1.SaveAs(f'{outputDir}/{histrdf}/{histcoll}/{histname}.png')

    if plotkey == 'autoCompPlot':
        if len(legList) != len(histList):
            raise RuntimeError("Legend list length unequal with length of histogram list!")

        histlist_copied = [hist.Clone() for hist in histList]

        histMetaName = histList[0].GetName()
        histrdf, histcoll, histname = histname_to_folder(histMetaName)
        drawableObjectList = makeGoodSinglePlotAuto(histlist_copied[0])

        rebinnedhistlist = autoModifyHists(histlist_copied, drawableObjectList[0])

        c1 = ROOT.TCanvas('c1', 'c1', 1000, 600)
        leg = ROOT.TLegend(0.0, 0.0, 1.0, 1.0)
        pad_plot = ROOT.TPad("pad_plot", "Left Pad", 0.0, 0.0, 0.8, 1.0)
        pad_leg = ROOT.TPad("pad_leg", "Right Pad", 0.8, 0.9-0.1*len(legList), 1.0, 0.9)
        pad_plot.Draw()
        pad_leg.Draw()

        pad_plot.cd()
        ROOT.gStyle.SetOptStat(0)

        logy, _ = auto_ylog_decision(drawableObjectList[0])
        pad_plot.SetLogy(logy)
        
        # if logy:
        #     drawableObjectList[0].SetMinimum(0.9)

        colours = generateColorPalette(len(legList))
        drawableObjectList[0].SetLineWidth(0)
        drawableObjectList[0].SetFillColorAlpha(colours[0], 0.25)
        name = drawableObjectList[0].GetName()
        scalefactor = 1.0/normlist[0] if not(name.endswith('cumulative')) else 1.0
        drawableObjectList[0].Scale(scalefactor)

        y_max_0 = max([drawableObjectList[0].GetBinContent(i) for i in range(drawableObjectList[0].GetNbinsX()+2)]) 
        y_max_list = [max([h.GetBinContent(i) for i in range(h.GetNbinsX())])/normlist[i+1] for i,h in enumerate(rebinnedhistlist)]
        y_max_list.append(y_max_0)
        y_max = max(y_max_list)
        # y_max = (15 if logy else 1.5)*y_max if y_max < 0.25 else (2 if logy else 1.1)
        drawableObjectList[0].SetMaximum(15*y_max if logy else 1.1*y_max)
        
        drawableObjectList[0].Draw('HIST E1')
        leg.AddEntry(drawableObjectList[0], legList[0], 'f')

        for i, hist in enumerate(rebinnedhistlist):
            hist.SetLineColorAlpha(colours[i+1], 0.35)
            name = hist.GetName()
            scalefactor = 1.0 / normlist[i+1] if not(name.endswith('cumulative')) and hist.Integral()!=0 else 1.0
            hist.Scale(scalefactor)
            hist.Draw('HIST E1 SAME')
            leg.AddEntry(hist, legList[i+1], 'lep')

        for i in range(1, len(drawableObjectList)):
            drawableObjectList[i].Draw()

        for canvconds in cc['canvas']:
            res = canvconds[0](histname)
            if res:
                canvconds[1](pad_plot)
        
        pad_leg.cd()
        leg.SetBorderSize(0)
        leg.Draw()

        if not os.path.exists(f'{outputDir}/{histrdf}/{histcoll}'):
            os.makedirs(f'{outputDir}/{histrdf}/{histcoll}')
        c1.SaveAs(f'{outputDir}/{histrdf}/{histcoll}/{histname}.png')
