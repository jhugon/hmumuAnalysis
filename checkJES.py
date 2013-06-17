import ROOT as root

def getHist(tree,var,cuts,name,nbins=50,minX=110,maxX=160):
   drawStr = var+" >> "+name+"("+str(nbins)+","+str(minX)+","+str(maxX)+")"
   tree.Draw(drawStr,cuts)
   tmp = root.gDirectory.Get(name)
   if type(tmp) != root.TH1F:
     print("Warning: loading histogram: '%s' from var '%s': Object type is not TH1F!!" % (name,var))
     print("  Draw string: '%s'" % drawStr)
     print("  cut string: '%s'" % cuts)
   #tmp.UseCurrentStyle()
   tmp.Sumw2()
   tmp.SetTitle("")
   return tmp

#f = root.TFile("vbfHmumu125_8TeV.root")
f = root.TFile("ggHmumu125_8TeV.root")

tree = f.Get('outtree')
canvas = root.TCanvas("c")

hists = []
for hName in ["dijetMass","deltaEtaJets","ptMiss","jetLead_pt","jetSub_eta","nJets"]:
  histsTmp = []
  nbins = 20
  minX = 0.0
  maxX = 500.
  maxY = 0.
  if hName == "nJets":
    nbins = 11
    minX = 0.0
    maxX = 10.
  elif "deltaEta" in hName:
    minX = 0.
    maxX = 10.
  elif "Eta" in hName or "eta" in hName:
    minX = -5.
    maxX = 5.
  for i,c in zip(["","_JESUp","_JESDown","_JERUp","_JERDown"],[1,root.kRed,root.kBlue,root.kRed,root.kBlue]):
    cuts = "nJets"+i+">=2"
    #cuts = "jetLead_pt"+i+">=30. && jetSub_pt"+i+">=30."
    if hName == "nJets":
      cuts = ""
    h = getHist(tree,hName+i,cuts,hName+i,nbins,minX,maxX)
    h.SetLineColor(c)
    h.SetMarkerColor(c)
    histsTmp.append(h)
    maxY = max(maxY,h.GetMaximum())
  histsTmp[0].GetYaxis().SetRangeUser(0.,maxY*1.05)
  histsTmp[0].Draw()
  histsTmp[1].Draw('same')
  histsTmp[2].Draw('same')
  canvas.SaveAs(hName+"JES.png")
  histsTmp[0].Draw()
  histsTmp[3].Draw('same')
  histsTmp[4].Draw('same')
  canvas.SaveAs(hName+"JER.png")
