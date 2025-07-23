import math as m
import numpy as np
import re
import os

import ROOT
from ROOT import TH1

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


canvas_customisation_conds = [
    (lambda histname: 'Iso' in histname, lambda canvas: canvas.SetLogx())]


def makePlot_isolation(histo : TH1):
    # Generate 100 xbins upto precision 2 in log-scale
    xlow = 0.001
    xhigh = 100
    nbins = 200
    step = (m.log10(xhigh) - m.log10(xlow))/nbins
    bins = [np.round(10**i, 3) for i in np.arange(m.log10(xlow), m.log10(xhigh), step)]
    bins = list(np.sort(list(set(bins))))
    bins.append(100)
    bins_array = np.array(bins, dtype='float64')
    histo.SetBinContent(2, histo.GetBinContent(1))
    histo.SetBinContent(1, 0)
    histo_rebinned = histo.Rebin(len(bins)-1, histo.GetName()+'_rebinned', bins_array)
    for bin in range(1, histo_rebinned.GetNbinsX()+1):
        binc = histo_rebinned.GetBinContent(bin)
        bw = np.round(histo_rebinned.GetBinWidth(bin), 3)
        histo_rebinned.SetBinContent(bin, binc/bw)
    histo_rebinned.GetYaxis().SetTitle('Events / bin width')
    return histo_rebinned


def defineXrangeAuto(histo : TH1, includeUnderflow : bool = True, includeOverflow : bool = True):
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


def makeBinLabel(histo : TH1, binNumber : int, label : str):
    x = histo.GetXaxis().GetBinLowEdge(binNumber)
    y = histo.GetBinContent(binNumber)
    labelobj = ROOT.TText(x, y, label)
    return labelobj


def makeGoodSinglePlotAuto(histo : TH1):
    if 'iso' in histo.GetTitle():
        histo = makePlot_isolation(histo)
    
    objectList = [histo]

    # The optimal x-range is defined as the first and the last bin
    # with non-zero entries including underflow and overflow
    binDecisions = defineXrangeAuto(histo)

    histo.SetLineWidth(4)
    if(binDecisions['uFlow']):
        ulabel = makeBinLabel(histo, 0, 'UFLOW')
        ulabel.SetTextAngle(90)
        ulabel.SetTextAlign(13)
        objectList.append(ulabel)
    if(binDecisions['oFlow']):
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
    tolerance = 5
    binentries = [hist.GetBinContent(i) for i in range(hist.GetNbinsX())]
    nonzerobinentries = [binentry for binentry in binentries if binentry != 0]
    binentrymax = max(nonzerobinentries)
    binentrymin = min(nonzerobinentries)
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
        hist_rebinned.SetLineWidth(4)
        hist_rebinned.GetXaxis().SetTitle(hist_rebinned.GetTitle())
        hist_rebinned.SetTitle('')
        newhistlist.append(hist_rebinned)

    return newhistlist


def generateColorPalette(ncolours):

    ROOT.gStyle.SetPalette(ROOT.kFruitPunch)

    colours = [ROOT.gStyle.GetColorPalette(int((i+1) * 255/ncolours)) for i in range(ncolours)]
    return colours


def makePngPlot(histList, outputDir:str, plotkey:str, legList=[]):

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

        for canvconds in canvas_customisation_conds:
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
        drawableObjectList[0].SetLineColor(colours[0])
        drawableObjectList[0].DrawNormalized('HIST E1')
        leg.AddEntry(drawableObjectList[0], legList[0], 'lep')

        for i, hist in enumerate(rebinnedhistlist):
            hist.SetLineColor(colours[i+1])
            hist.DrawNormalized('HIST E1 SAME')
            leg.AddEntry(hist, legList[i+1], 'lep')

        for i in range(1, len(drawableObjectList)):
            drawableObjectList[i].Draw()

        for canvconds in canvas_customisation_conds:
            res = canvconds[0](histname)
            if res:
                canvconds[1](pad_plot)

        pad_leg.cd()
        leg.SetBorderSize(0)
        leg.Draw()

        if not os.path.exists(f'{outputDir}/{histrdf}/{histcoll}'):
            os.makedirs(f'{outputDir}/{histrdf}/{histcoll}')
        c1.SaveAs(f'{outputDir}/{histrdf}/{histcoll}/{histname}.png')
