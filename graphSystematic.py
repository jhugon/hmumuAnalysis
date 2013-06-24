#!/usr/bin/env python

import os
import sys
import re

import ROOT as root
from ROOT import gStyle as gStyle
root.gErrorIgnoreLevel = root.kWarning
root.gROOT.SetBatch(True)

def setStyle():
  gStyle.SetCanvasColor(0)
  gStyle.SetCanvasBorderSize(10)
  gStyle.SetCanvasBorderMode(0)
  gStyle.SetCanvasDefH(700)
  gStyle.SetCanvasDefW(700)

  gStyle.SetPadColor       (0)
  gStyle.SetPadBorderSize  (10)
  gStyle.SetPadBorderMode  (0)
  gStyle.SetPadBottomMargin(0.13)
  gStyle.SetPadTopMargin   (0.08)
  gStyle.SetPadLeftMargin  (0.15)
  gStyle.SetPadRightMargin (0.05)
  gStyle.SetPadGridX       (0)
  gStyle.SetPadGridY       (0)
  gStyle.SetPadTickX       (1)
  gStyle.SetPadTickY       (1)

  gStyle.SetFrameFillStyle ( 0)
  gStyle.SetFrameFillColor ( 0)
  gStyle.SetFrameLineColor ( 1)
  gStyle.SetFrameLineStyle ( 0)
  gStyle.SetFrameLineWidth ( 1)
  gStyle.SetFrameBorderSize(10)
  gStyle.SetFrameBorderMode( 0)

  gStyle.SetNdivisions(505)

  gStyle.SetLineWidth(2)
  gStyle.SetHistLineWidth(2)
  gStyle.SetFrameLineWidth(2)
  gStyle.SetLegendFillColor(root.kWhite)
  gStyle.SetLegendFont(42)
  gStyle.SetMarkerSize(1.2)
  gStyle.SetMarkerStyle(20)
 
  gStyle.SetLabelSize(0.040,"X")
  gStyle.SetLabelSize(0.040,"Y")

  gStyle.SetLabelOffset(0.010,"X")
  gStyle.SetLabelOffset(0.010,"Y")
 
  gStyle.SetLabelFont(42,"X")
  gStyle.SetLabelFont(42,"Y")
 
  gStyle.SetTitleBorderSize(0)
  gStyle.SetTitleFont(42)
  gStyle.SetTitleFont(42,"X")
  gStyle.SetTitleFont(42,"Y")

  gStyle.SetTitleSize(0.045,"X")
  gStyle.SetTitleSize(0.045,"Y")
 
  gStyle.SetTitleOffset(1.4,"X")
  gStyle.SetTitleOffset(1.4,"Y")
 
  gStyle.SetTextSize(0.055)
  gStyle.SetTextFont(42)
 
  gStyle.SetOptStat(0)
  
setStyle()

def setHistTitles(hist,xlabel,ylabel):
    hist.GetXaxis().SetTitle(xlabel)
    hist.GetYaxis().SetTitle(ylabel)

def getUncertainty(nomFilename,sysFileNames):
  effFile = open(nomFilename, 'r')
  sysFiles = [open(i, 'r') for i in sysFileNames]
  
  matchString = r"([\w]+)[\s]+([-0-9.eE]+)[\s]+([-0-9.eE]+)[\s]*"

  uncertainties = {}
  statErrors = {}
  
  while True:
    try:
      nomLine = effFile.next()
      sysLines = [i.next() for i in sysFiles]
      if not "Jet" in nomLine: # to only do Baseline++ cats
        continue
      nomMatch = re.match(matchString,nomLine)
      if not nomMatch:
        print("Error: Nom line didn't match: %s" % nomLine)
        continue
      nomName = nomMatch.group(1)
      nomEff = float(nomMatch.group(2))
      statErr = float(nomMatch.group(3))
      relStatErr = statErr/nomEff
      errEffs = []
      correlation = None
      for i in sysLines:
        errMatch = re.match(matchString,i)
        if not errMatch:
          print("Error: Err line didn't match: %s" % i)
          errEffs.append(0.)
          continue
        tmpErrName = errMatch.group(1)
        tmpErrEff = float(errMatch.group(2))
        if tmpErrName != tmpErrName:
          print("Error: Err cat name didn't match nom: %s != %s" % (tmpErrName, nomName))
          errEff.append(0.)
          continue
        errEffs.append(abs(tmpErrEff-nomEff))
        if correlation == None:
          if tmpErrEff>=nomEff:
            correlation = ""
          else:
            correlation = "-"
          
      relErr = max(errEffs)/nomEff
      #print "%s  %s%.4f" % (nomName, correlation,relErr)
      #print "%s  %s%.4f  %.4f" % (nomName, correlation, relErr,relStatErr)
      relErr = abs(relErr)
      if correlation == "-":
        relErr *= -1
      uncertainties[nomName] = relErr
      statErrors[nomName] = relStatErr
    except StopIteration:
      break
  
  # cleanup
  for i in sysFiles:
    i.close()
  effFile.close()
  return uncertainties, statErrors

energies = ["8TeV"]
masses = ["125"]
energies = ["8TeV","7TeV"]
masses = ["115","125","135","150"]
errorSets = {
  "Scale Variations":["0p5X","2X"],
  "UE Variations":["D6T","P0","ProPT0","ProQ20","UEOFF"]
}
categories = [
  "Jets01PassPtG10",
  "Jets01FailPtG10",
  "Jet2CutsVBFPass",
  "Jet2CutsGFPass",
  "Jet2CutsFailVBFGF"
]
labels = [
  "Non-VBF Tight",
  "Non-VBF Loose",
  "VBF-Presel VBF Tight",
  "VBF-Presel GF Tight",
  "VBF-Presel Loose"
]
colors = [
  root.kBlue,
  root.kRed,
  root.kGreen,
  root.kMagenta,
  root.kCyan
]
shifts = [0.2*(-len(categories)/2.+i) for i in range(len(categories))]

canvas = root.TCanvas("c1")

for errorSet in errorSets:
  dataMasses = {}
  dataUnc = {}
  dataStat = {}
  for energy in energies:
    tmpMassesList = []
    tmpUncList = []
    tmpStatList = []
    for mass in sorted(masses,key=float):
      nomFn =  "ggHmumu%s_%s.root.txt" % (mass,energy)
      if not os.path.exists(nomFn):
        print("Nominal File doesn't exist: %s" % nomFn)
        continue
      errFns =  ["ggHmmu%s%s_%s.root.txt" % (energy,mass,i) for i in errorSets[errorSet]]
      sysNotExist = False
      for fn in errFns:
        if not os.path.exists(fn):
          print("Systematic File doesn't exist: %s" % fn)
          sysNotExist = True
      if sysNotExist:
        continue
      tmpUnc, tmpStat = getUncertainty(nomFn, errFns)
      tmpMassesList.append(float(mass))
      tmpUncList.append(tmpUnc)
      tmpStatList.append(tmpStat)
    dataMasses[energy] = tmpMassesList
    dataUnc[energy] = tmpUncList
    dataStat[energy] = tmpStatList

  graphs = []
  axisList = []
  for energy in energies:
    axis = root.TH2F("axis"+errorSet+energy,"",1,110,160,1,0.0,50)
    axis.GetXaxis().SetTitle("m_{H} [GeV/c^{2}]")
    axis.GetYaxis().SetTitle("Uncertainty on Signal Yield [%]")
    axis.Draw()
    axisList.append(axis)
    leg = root.TLegend(0.58,0.70,0.9,0.9)
    #leg = root.TLegend(0.45,0.60,0.95,0.9)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    for cat,col,shift,label in zip(categories,colors,shifts,labels):
      graph = root.TGraphErrors()
      graph.SetLineColor(col)
      graph.SetMarkerColor(col)
      for iMass in range(len(dataMasses[energy])):
        x = dataMasses[energy][iMass]
        y = dataUnc[energy][iMass][cat]
        yErr = dataStat[energy][iMass][cat]*100.
        y = abs(y)*100.
        #print("x = %s, y = %s" % (x,y))
        graph.SetPoint(iMass,x+shift,y)
        graph.SetPointError(iMass,0.,yErr)
      graph.Draw("PEZ")
      #graph.Print()
      graphs.append(graph)
      leg.AddEntry(graph,label,"P")
    leg.Draw()
    errorSetSaveName = errorSet.replace(" ","") 
    canvas.SaveAs(errorSetSaveName+"_"+energy+".png")
