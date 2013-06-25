#!/usr/bin/env python

import os
import sys
import re

import ROOT as root
root.gROOT.SetBatch(True)
root.gErrorIgnoreLevel = root.kWarning
from ROOT import gStyle as gStyle

THRESHOLD=0.01
PRELIMINARYSTRING="CMS Preliminary"

def absMax(x,y):
  absx = abs(x)
  absy = abs(y)
  if absx > absy:
    return x
  else:
    return y

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

def getEfficiencies(filename):
  effFile = open(filename, 'r')
  matchString = r"([\w]+)[\s]+([-0-9.eE]+)[\s]+([-0-9.eE]+)[\s]*"
  effs = {}
  for line in effFile:
    nomMatch = re.match(matchString,line)
    if not nomMatch:
      print("Error: Eff line didn't match: %s" % line)
      continue
    nomName = nomMatch.group(1)
    nomEff = float(nomMatch.group(2))
    statErr = float(nomMatch.group(3))
    effs[nomName] = {"eff":float(nomEff),"err":float(statErr)}
  effFile.close()
  return effs

#datasets = ["ggH"]
#energies = ["8TeV"]
#masses = ["125"]
datasets = ["ggH","vbfH"]
energies = ["8TeV","7TeV"]
masses = ["115","125","135","150"]
errorSets = {
  "Scale Variations":["0p5X","2X"],
  "UE Variations":["D6T","P0","ProPT0","ProQ20","UEOFF"]
}
errorSetKeyName = {
  "Scale Variations":"scale",
  "UE Variations":"ue"
}
categories = [
  "Jets01PassPtG10",
  "Jets01FailPtG10",
  "Jet2CutsVBFPass",
  "Jet2CutsGFPass",
  "Jet2CutsFailVBFGF"
]
labels = [
  "Non-VBF Presel. Tight",
  "Non-VBF Presel. Loose",
  "VBF-Presel. VBF Tight",
  "VBF-Presel. GF Tight",
  "VBF-Presel. Loose"
]
colors = [
  root.kBlue,
  root.kRed,
  root.kGreen,
  root.kMagenta,
  root.kCyan,
  root.kOrange+3,
  root.kOrange-3,
  root.kGray+2
]
shifts = [0.2*(-len(categories)/2.+i) for i in range(len(categories))]

canvas = root.TCanvas("c1")
tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.04)
tlatex.SetTextAlign(12)

for errorSet in errorSets:
  dataMasses = {}
  dataUnc = {}
  dataStat = {}
  for energy in energies:
    dataMasses[energy] = {}
    dataUnc[energy] = {}
    dataStat[energy] = {}
    for ds in datasets:
      tmpMassesList = []
      tmpUncList = []
      tmpStatList = []
      for mass in sorted(masses,key=float):
        nomFn =  "%smumu%s_%s.root.txt" % (ds,mass,energy)
        if not os.path.exists(nomFn):
          print("Nominal File doesn't exist: %s" % nomFn)
          continue
        errFns =  ["%smmu%s%s_%s.root.txt" % (ds,energy,mass,i) for i in errorSets[errorSet]]
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
      dataMasses[energy][ds] = tmpMassesList
      dataUnc[energy][ds] = tmpUncList
      dataStat[energy][ds] = tmpStatList
  
  graphs = []
  axisList = []
  for energy in energies:
    for ds in datasets:
      axis = root.TH2F("axis"+errorSet+energy+ds,"",1,110,160,1,0.0,50)
      axis.GetXaxis().SetTitle("m_{H} [GeV/c^{2}]")
      axis.GetYaxis().SetTitle("Uncertainty on Signal Yield [%]")
      axis.Draw()
      axisList.append(axis)
      leg = root.TLegend(0.58,0.70,0.9,0.9)
      #leg = root.TLegend(0.45,0.60,0.95,0.9)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      for cat,col,shift,label in zip(categories,colors[:len(categories)],shifts,labels):
        graph = root.TGraphErrors()
        graph.SetLineColor(col)
        graph.SetMarkerColor(col)
        for iMass in range(len(dataMasses[energy][ds])):
          x = dataMasses[energy][ds][iMass]
          y = dataUnc[energy][ds][iMass][cat]
          yErr = dataStat[energy][ds][iMass][cat]*100.
          y = abs(y)*100.
          #print("x = %s, y = %s" % (x,y))
          graph.SetPoint(iMass,x+shift,y)
          graph.SetPointError(iMass,0.,yErr)
        graph.Draw("PEZ")
        #graph.Print()
        graphs.append(graph)
        leg.AddEntry(graph,label,"P")
      leg.Draw()
      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
      tlatex.SetTextAlign(32)
      tlatex.DrawLatex(1.02-gStyle.GetPadRightMargin(),0.96,errorSet)
      tlatex.SetTextAlign(12)
      dsLabel = "GF"
      if ds=="vbfH":
        dsLabel = "VBF"
      captionStr = r"%s %s" % (dsLabel,energy.replace("TeV"," TeV"))
      tlatex.DrawLatex(0.04+gStyle.GetPadLeftMargin(),0.88,captionStr)
      errorSetSaveName = errorSet.replace(" ","") 
      canvas.SaveAs(errorSetSaveName+"_"+ds+energy+".png")
  
  # Now for a table
  tableStr = ""
  #extraCols = len(datasets)
  #tableStr += r"\begin{tabular}{|c|l|"+"c|"*extraCols+r"} \hline" + "\n"
  #tableStr += r"$\sqrt{s}$ & Category & "
  #for ds in datasets:
  #  dsLabel = "GF"
  #  if ds=="vbfH":
  #    dsLabel = "VBF"
  #  tableStr += r"%s &" % dsLabel
  #tableStr = tableStr[:-1] + r"\\ \hline \hline" + "\n"
  #for energy in energies:
  #  tableStr += "\multirow{"+str(len(categories))+"}{*}"
  #  tableStr += "{%s} \n" % energy.replace("TeV"," TeV")
  #  for cat,label in zip(categories,labels):
  #    tableStr += " & "+label + " & "
  #    for ds in datasets:
  #      yMax = 0.
  #      for iMass in range(len(dataMasses[energy][ds])):
  #        x = dataMasses[energy][ds][iMass]
  #        y = dataUnc[energy][ds][iMass][cat]
  #        yMax = absMax(y,yMax)
  #      tableStr +=  "%.2f%% &" % (yMax*100.)
  #    tableStr = tableStr[:-1] + r" \\ "+ "\n"
  #  tableStr = tableStr[:-1] + r" \hline "+ "\n"
  #tableStr += r"\end{tabular}" + "\n"

  extraCols = len(datasets)*len(energies)
  tableStr += r"\begin{tabular}{|l|"+"c|"*extraCols+r"} \hline" + "\n"
  tableStr += r"Category &"
  for ds in datasets:
    dsLabel = "GF"
    if ds=="vbfH":
      dsLabel = "VBF"
    for energy in energies:
      tableStr += r" %s %s &" % (dsLabel,energy.replace("TeV"," TeV"))
  tableStr = tableStr[:-1] + r"\\ \hline \hline" + "\n"
  for cat,label in zip(categories,labels):
    tableStr += label + " &"
    for ds in datasets:
      for energy in energies:
        yMax = 0.
        for iMass in range(len(dataMasses[energy][ds])):
          x = dataMasses[energy][ds][iMass]
          y = dataUnc[energy][ds][iMass][cat]
          yMax = absMax(y,yMax)
        tableStr +=  " %.2f%% &" % (yMax*100.)
    tableStr = tableStr[:-1] + r" \\ \hline "+ "\n"
  tableStr += r"\end{tabular}" + "\n"

  tableStr = tableStr.replace(r"%",r"\%")
  print
  print errorSet
  print
  print tableStr
  print

  print("    self."+errorSetKeyName[errorSet]+" = {")
  for dsName in datasets:
    dsLabel = "gg"
    if ds=="vbfH":
      dsLabel = "vbf"
    print "      '"+dsLabel+"' : {"
    for energy in energies:
      print "        '"+energy+"' : {"
      for key in categories:
        maxErr = 0.
        for iMass in range(len(dataMasses[energy][ds])):
          x = dataMasses[energy][dsName][iMass]
          y = dataUnc[energy][dsName][iMass][key]
          maxErr = absMax(y,maxErr)
        if abs(maxErr) < THRESHOLD:
          print "          '"+key+"' : None,"
        else:
          if maxErr>0.:
            maxErr += 1.
          else:
            maxErr -= 1.
          print "          '"+key+"' : %.4f," % (maxErr)
      print("          },")
    print("        },")
  # for wH and zH
  for dsName in datasets:
    dsPrint = "w"
    if "vbf" in dsName:
      dsPrint = "z"
    print "      '"+dsPrint+"' : {"
    for energy in energies:
      print "        '"+energy+"' : {"
      for key in categories:
          print "          '"+key+"' : None,"
      print("          },")
    print("        },")
  print("      }")
  print
  
# Now UE specific plot
ueMasses = masses
ueMasses.remove("135")
for ds in datasets:
  data = {}
  for energy in energies:
    data[energy] = {}
    for mass in ueMasses:
      data[energy][mass] = {}
      # Also do nominal
      nomFn =  "%smumu%s_%s.root.txt" % (ds,mass,energy)
      if not os.path.exists(nomFn):
        #print("Nominal File doesn't exist: %s" % nomFn)
        continue
      nomEffs = getEfficiencies(nomFn)
      data[energy][mass]["nom"] = nomEffs
      # Now do variations
      for error in errorSets["UE Variations"]:
        errFn =  "%smmu%s%s_%s.root.txt" % (ds,energy,mass,error)
        if not os.path.exists(errFn):
          continue
        errEffs = getEfficiencies(errFn)
        data[energy][mass][error] = errEffs
  for cat,label in zip(categories,labels):
    graphs = []
    leg = root.TLegend(0.75,0.66,0.9,0.9)
    #leg = root.TLegend(0.45,0.60,0.95,0.9)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    yMax = 0.
    yMin = 1.
    for energy in energies:
      for error,color in zip(["nom"]+errorSets["UE Variations"],[1]+colors[:len(errorSets["UE Variations"])]):
        graph = root.TGraphErrors()
        iPoint = 0
        for mass in masses:
          if not data[energy].has_key(mass):
            continue
          if not data[energy][mass].has_key(error):
            continue
          tmpMass = float(mass)
          tmpEff = data[energy][mass][error][cat]["eff"]
          tmpErr = data[energy][mass][error][cat]["err"]
          tmpEff *= 100.
          tmpErr *= 100.
          graph.SetPoint(iPoint,tmpMass,tmpEff)
          graph.SetPointError(iPoint,0.,tmpErr)
          yMin = min(yMin,tmpEff-tmpErr)
          yMax = max(yMax,tmpEff+tmpErr)
          iPoint += 1
        graph.SetMarkerColor(color)
        graph.SetLineColor(color)
        #print "energy: %s mass: %.1f error: %s cat: %s" % (energy,tmpMass,error,cat)
        #print "  %.2f +/- %.2e" % (tmpEff,tmpErr)
        graph.Draw("ap")
        if energy=="8TeV":
          if error == "nom":
            error = "Z2 / Z2*"
          leg.AddEntry(graph,error,"lep")
        else:
          graph.SetLineWidth(graph.GetLineWidth()/2)
          graph.SetLineStyle(2)
          graph.SetMarkerStyle(24)
        graphs.append(graph)
    canvas.Clear()
    #axis = root.TH2F("axisTunes"+cat+ds+energy,"",1,110,160,1,0.0,0.5)
    axis = root.TH2F("axisTunes"+cat+ds+energy,"",1,110,160,1,0.,yMax*2.)
    axis.GetXaxis().SetTitle("m_{H} [GeV/c^{2}]")
    axis.GetYaxis().SetTitle("Signal Efficiency [%]")
    axis.Draw()
    for graph in graphs:
      graph.Draw("PL")
    leg.Draw()
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.02-gStyle.GetPadRightMargin(),0.96,label)
    tlatex.SetTextAlign(12)
    dsLabel = "GF"
    if ds=="vbfH":
      dsLabel = "VBF"
    captionStr = r"%s" % (dsLabel)
    tlatex.DrawLatex(0.04+gStyle.GetPadLeftMargin(),0.88,captionStr)
    errorSetSaveName = errorSet.replace(" ","") 
    canvas.SaveAs("TuneEffComp_"+cat+"_"+ds+".png")
