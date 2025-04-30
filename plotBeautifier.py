from datetime import datetime
import numpy as np
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
            beginBin = i
            break

    # Loop over the bin content and find the last bin with data in it
    endBin = histo.GetNbinsX()+1 if includeOverflow else histo.GetNbinsX()
    for i in range(endBin, beginBin-1, -1):
        if histo.GetBinContent(i) != 0:
            endBin = i
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

    binDecisions = defineXrangeAuto(histo)

    histo.SetTitle('')
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

    return objectList


def makePngPlot(histObj:TH1, outputDir:str, plotkey:str):
    histName = histObj.GetName()
    histCopy = histObj.Clone()

    if plotkey == 'autoSinglePlot':
        drawableObjectList = makeGoodSinglePlotAuto(histCopy)

    c1 = ROOT.TCanvas('c1', 'c1', 800, 600)
    ROOT.gStyle.SetOptStat(0)
    drawableObjectList[0].Draw()

    for i in range(1, len(drawableObjectList)):
        drawableObjectList[i].Draw()

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    c1.SaveAs(f'{outputDir}/{histName}.png')
    