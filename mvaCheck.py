#!/usr/bin/python

import os.path
import re
import glob
import ROOT as root

root.gROOT.SetBatch(True)
tools = root.TMVA.Tools.Instance();
root.gStyle.SetPaintTextFormat("3g")

def printSeparations(infilename):
  separations = {}
  fileNameMatch = re.match("TMVA_(.+)_(.+).root",infilename)
  infile = root.TFile(infilename)
  for fileKey in infile.GetListOfKeys():
    fileKeyName = fileKey.GetName()
    #match = re.match("Method_(.+)",fileKeyName)
    match = re.match("Method_BDT",fileKeyName)
    if match:
      folder1 = fileKey.ReadObj()
      for folder1Key in folder1.GetListOfKeys():
        folder2 = folder1Key.ReadObj()
        filePathName = folder1.GetName()+"/"+folder2.GetName()
        varNames = set()
        for folder2Key in folder2.GetListOfKeys():
          keyName =  folder2Key.GetName()
          matchS = re.match(r"(.+)__Signal",keyName)
          #matchB = re.match(r"(.+)__Background",keyName)
          if matchS:
            varNames.add(matchS.group(1))
        for varName in varNames:
          sHist = folder2.Get(varName+"__Signal")
          bHist = folder2.Get(varName+"__Background")
          sep = tools.GetSeparation(sHist,bHist)
          separations[varName] = sep
        sHist = folder2.Get("MVA_BDT_S")
        bHist = folder2.Get("MVA_BDT_B")
        sep = tools.GetSeparation(sHist,bHist)
        separations["BDT"] = sep
        
  fileTitle = "Unknown file!"
  if fileNameMatch:
    if fileNameMatch.group(1) == "vbf":
      fileTitle = "VBF"
    elif fileNameMatch.group(1) == "inclusive":
      fileTitle = "Non-VBF"
    energyStr = fileNameMatch.group(2)
    fileTitle += energyStr.replace("TeV"," TeV")
  print(fileTitle)
  print("-"*40)
  for i in reversed(sorted(separations.keys(),key=lambda x: separations[x])):
    print("%-20s: %g" % (i,separations[i]))
  print
  canvas = root.TCanvas()
  fileTitle = fileTitle.replace(" ","")
  fileTitle = fileTitle.replace("-","")
  for i in ["S","B"]:
    objName = "CorrelationMatrix"+i
    corMat = infile.Get(objName)
    #corMat.Scale(0.01)
    corMat.Draw("coltext")
    canvas.SaveAs(fileTitle+"_"+objName+".png")
    canvas.SaveAs(fileTitle+"_"+objName+".pdf")
    #canvas.SaveAs(fileTitle+"_"+objName+".root")
    #canvas.SaveAs(fileTitle+"_"+objName+".eps")

if __name__ == "__main__":
  for i in glob.glob("TMVA*.root"):
    printSeparations(i)
