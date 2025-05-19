import re
import os

import ROOT
from ROOT import TH1

def defineXrangeAuto(histo : TH1, includeUnderflow : bool = True, includeOverflow : bool = True):
    # Decide if underflow and overflow bins should be included in the plot
    includeUnderflow = histo.GetBinContent(0) > 0.01*histo.GetBinContent(1)
    includeOverflow = histo.GetBinContent(histo.GetNbinsX()+1) > 0.01*histo.GetBinContent(histo.GetNbinsX())
    
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


def makePngPlot(histObj:TH1, outputDir:str, plotkey:str):
    histMetaName = histObj.GetName()
    histCopy = histObj.Clone()

    # Get the folder and histogram details from the histogram name
    # Assuming name has been saved as DATAFRAME_COLLECTION_<OBJSELECTION_>VAR
    # The OBJSELECTION could be optional
    histMetaNameSplit = re.fullmatch('([A-Za-z0-9]+)_([A-Za-z0-9]+)_([A-Za-z0-9_]+)', histMetaName)
    histrdf = histMetaNameSplit.group(1)
    histcoll = histMetaNameSplit.group(2)
    histname = histMetaNameSplit.group(3)

    if plotkey == 'autoSinglePlot':
        drawableObjectList = makeGoodSinglePlotAuto(histCopy)

    c1 = ROOT.TCanvas('c1', 'c1', 800, 600)
    ROOT.gStyle.SetOptStat(0)

    # If the minimum non-zero bin entry is less than tolerance percentage
    # of the maximum bin entry, plot in log scale
    tolerance = 1
    binentries = [drawableObjectList[0].GetBinContent(i) for i in range(drawableObjectList[0].GetNbinsX())]
    nonzerobinentries = [binentry for binentry in binentries if binentry != 0]
    binentrymax = max(nonzerobinentries)
    binentrymin = min(nonzerobinentries)
    if binentrymin < tolerance * 0.01 * binentrymax:
        c1.SetLogy()

    drawableObjectList[0].Draw()

    for i in range(1, len(drawableObjectList)):
        drawableObjectList[i].Draw()

    if not os.path.exists(f'{outputDir}/{histrdf}/{histcoll}'):
        os.makedirs(f'{outputDir}/{histrdf}/{histcoll}')
    c1.SaveAs(f'{outputDir}/{histrdf}/{histcoll}/{histname}.png')
    