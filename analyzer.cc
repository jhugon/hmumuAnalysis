#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "DataFormats.h"
#include "helpers.h"
#include "mva.h"
#include "LumiReweightingStandAlone.h"

#include "boost/program_options.hpp"
#include "boost/regex.hpp"

#include <limits.h>

#define JETPUID
//#define PUREWEIGHT

using namespace std;
using namespace boost;

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  const char* optionIntro = "H->MuMu Analyzer\n\nUsage: ./analyzer [--help] [--train] [--maxEvents N] <outputFileName.root> <inputFileName.root> [<inputFileName2.root>...]\n\nAllowed Options";
  program_options::options_description optionDesc(optionIntro);
  optionDesc.add_options()
      ("help,h", "Produce Help Message")
      ("trainingTree,t",program_options::value<string>(), "Create Training Tree File with filename")
      ("filenames",program_options::value<vector<string> >(), "Input & Output File Names, put output name first followed by all input file names")
      ("maxEvents,m",program_options::value<int>(), "Maximum Number of Events to Process")
  ;
  
  program_options::positional_options_description optionPos;
  optionPos.add("filenames",-1);
  
  program_options::variables_map optionMap;
  program_options::store(program_options::command_line_parser(argc, argv).options(optionDesc).positional(optionPos).run(), optionMap);
  program_options::notify(optionMap);    
  
  if (optionMap.count("help")) 
  {
      cout << optionDesc << "\n";
      return 1;
  }

  std::vector<std::string> filenames;
  vector<string>::const_iterator filename;
  std::string outputFileName;
  if (optionMap.count("filenames")>0)
  {
     filenames = optionMap["filenames"].as<vector<string> >();
     if(filenames.size()<2)
     {
       cout << "Error: Need both input file and output file names, exiting." << endl;
       return 1;
     }
     outputFileName = filenames[0];
     filenames.erase(filenames.begin());
  }
  else
  {
     cout << "Error: Input file name  and ouput file name arguments required, exiting." << endl;
     return 1;
  }

  int maxEvents = std::numeric_limits<int>::max();
  if (optionMap.count("maxEvents")) 
  {
      int tmp = optionMap["maxEvents"].as<int>();
      if (tmp > 0)
      {
        maxEvents = tmp;
      }
  }
  cout << "maxEvents = "<< maxEvents << "\n";

  /////////////////////////////
  //////////// Setup //////////
  /////////////////////////////

  float minMmm = 70.0;
  float maxMmm = 200.0;
  //float minMmm = 110.0;
  //float maxMmm = 150.0;

  float minBlind = 115;
  float maxBlind = 135;

  // Check to see if it is data
  bool isData = false;
  std::vector<std::string> dataWords;
  dataWords.push_back("Run2012");
  dataWords.push_back("Run2011");
  dataWords.push_back("SingleMu");
  dataWords.push_back("DoubleMu");
  std::vector<std::string>::const_iterator dataWord;
  for(dataWord = dataWords.begin(); dataWord != dataWords.end();dataWord++)
  {
    regex re(*dataWord);
    if(regex_search(filenames[0],re))
    {
        isData = true;
    }
  }

  if(isData)
    std::cout << "This is a Real Data Sample\n";
  else
    std::cout << "This is a MC Sample\n";

  ////////////
  
  TChain * tree = new TChain("tree");

  cout << "Input File Names: \n"; 
  for(filename = filenames.begin();filename != filenames.end();filename++)
  {
    cout<<"  "<< *filename << endl;
    tree->AddFile(filename->c_str());
  }

  cout << "Output File Name: " << outputFileName << endl;
  TFile * outFile = new TFile(outputFileName.c_str(),"RECREATE");

  std::string trainingTreeFileName = "";
  if (optionMap.count("trainingTree")) 
  {
      cout << "Training enabled" << "\n";
      trainingTreeFileName = optionMap["trainingTree"].as<string>();
      cout << "Training Tree File Name: " << trainingTreeFileName << "\n";
  }

  //////////////////////////
  // Tree Branches

  _MuonInfo reco1, reco2;

  tree->SetBranchAddress("reco1", &reco1);
  tree->SetBranchAddress("reco2", &reco2);

  float recoCandMass, recoCandPt, recoCandY, recoCandPhi;

  tree->SetBranchAddress("recoCandMass", &recoCandMass);
  tree->SetBranchAddress("recoCandPt"  , &recoCandPt );
  tree->SetBranchAddress("recoCandY"  , &recoCandY );
  tree->SetBranchAddress("recoCandPhi"  , &recoCandPhi );

  _PFJetInfo jets;
  tree->SetBranchAddress("pfJets",&jets);

#ifdef JETPUID
  float puJetFullDisc[10];
  float puJetSimpleDisc[10];
  float puJetCutDisc[10];

  tree->SetBranchAddress("puJetFullDisc",&puJetFullDisc);
  tree->SetBranchAddress("puJetSimpleDisc",&puJetSimpleDisc);
  tree->SetBranchAddress("puJetCutDisc",&puJetCutDisc);

  int puJetFullId[10];
  int puJetSimpleId[10];
  int puJetCutId[10];

  tree->SetBranchAddress("puJetFullId",&puJetFullId);
  tree->SetBranchAddress("puJetSimpleId",&puJetSimpleId);
  tree->SetBranchAddress("puJetCutId",&puJetCutId);
#endif

#ifdef PUREWEIGHT
  int nPU;
  tree->SetBranchAddress("nPU",&nPU);
#else
  int nPU=0;
#endif
  _VertexInfo vertexInfo;
  tree->SetBranchAddress("vertexInfo",&vertexInfo);
  _EventInfo eventInfo;
  tree->SetBranchAddress("eventInfo",&eventInfo);
  _MetInfo met;
  tree->SetBranchAddress("met",&met);

  //////////////////////////
  // Histograms
  std::map<std::string,TH1F*> histMap;
  std::map<std::string,TH2F*> histMap2D;
  std::map<std::string,TH1F*>::iterator histMapIter;
  std::map<std::string,TH2F*>::iterator histMap2DIter;

  TH1F* mDiMu = new TH1F("mDiMu","DiMuon Mass",1600,0,400);
  histMap.insert(make_pair("mDiMu",mDiMu));

  TH1F* mDiJet = new TH1F("mDiJet","DiJet Mass",500,0,2000);
  histMap.insert(make_pair("mDiJet",mDiJet));

  TH1F* ptDiMu = new TH1F("ptDiMu","DiMuon Pt",250,0,500);
  histMap.insert(make_pair("ptDiMu",ptDiMu));

  TH1F* yDiMu = new TH1F("yDiMu","DiMuon Rapidity",100,-4,4);
  histMap.insert(make_pair("yDiMu",yDiMu));

  TH2F* yVptDiMu = new TH2F("yVptDiMu","DiMuon Rapidity v. p_{T}",250,0,500,100,0,4);
  histMap2D.insert(make_pair("yVptDiMu",yVptDiMu));
  TH2F* ptVmDiMu = new TH2F("ptVmDiMu","DiMuon p_{T} v. Mass",1600,0,400,250,0,250);
  histMap2D.insert(make_pair("ptVmDiMu",ptVmDiMu));
  TH2F* yVmDiMu = new TH2F("yVmDiMu","DiMuon |y| v. Mass",1600,0,400,100,0,4);
  TH2F* phiVmDiMu = new TH2F("phiVmDiMu","DiMuon #phi v. Mass",1600,0,400,100,0,3.2);
  histMap2D.insert(make_pair("phiVmDiMu",phiVmDiMu));

  TH1F* ptMu1 = new TH1F("ptMu1","Leading Muon Pt",100,0,200);
  histMap.insert(make_pair("ptMu1",ptMu1));
  TH1F* ptMu2 = new TH1F("ptMu2","Sub-Leading Muon Pt",100,0,200);
  histMap.insert(make_pair("ptMu2",ptMu2));
  TH1F* ptJet1 = new TH1F("ptJet1","Leading Jet Pt",200,0,1000);
  histMap.insert(make_pair("ptJet1",ptJet1));
  TH1F* ptJet2 = new TH1F("ptJet2","Sub-Leading Jet Pt",200,0,1000);
  histMap.insert(make_pair("ptJet2",ptJet2));

  TH1F* etaMu1 = new TH1F("etaMu1","Leading Muon #eta",50,-2.5,2.5);
  histMap.insert(make_pair("etaMu1",etaMu1));
  TH1F* etaMu2 = new TH1F("etaMu2","Sub-Leading Muon #eta",50,-2.5,2.5);
  histMap.insert(make_pair("etaMu2",etaMu2));
  TH1F* etaJet1 = new TH1F("etaJet1","Leading Jet #eta",50,-5.0,5.0);
  histMap.insert(make_pair("etaJet1",etaJet1));
  TH1F* etaJet2 = new TH1F("etaJet2","Sub-Leading Jet #eta",50,-5.0,5.0);
  histMap.insert(make_pair("etaJet2",etaJet2));

  TH1F* deltaEtaJetsHist = new TH1F("deltaEtaJets","#Delta#eta Jets",50,0.0,10.0);
  histMap.insert(make_pair("deltaEtaJets",deltaEtaJetsHist));
  TH1F* deltaPhiJetsHist = new TH1F("deltaPhiJets","#Delta#phi Jets",50,0.0,3.2);
  histMap.insert(make_pair("deltaPhiJets",deltaPhiJetsHist));
  TH1F* deltaRJetsHist = new TH1F("deltaRJets","#Delta R Jets",50,0.0,10.0);
  histMap.insert(make_pair("deltaRJets",deltaRJetsHist));

  TH1F* deltaEtaMuonsHist = new TH1F("deltaEtaMuons","#Delta#eta Jets",50,0.0,10.0);
  histMap.insert(make_pair("deltaEtaMuons",deltaEtaMuonsHist));
  TH1F* deltaPhiMuonsHist = new TH1F("deltaPhiMuons","#Delta#phi Jets",50,0.0,3.2);
  histMap.insert(make_pair("deltaPhiMuons",deltaPhiMuonsHist));
  TH1F* deltaRMuonsHist = new TH1F("deltaRMuons","#Delta R Jets",50,0.0,10.0);
  histMap.insert(make_pair("deltaRMuons",deltaRMuonsHist));

  TH1F* countsHist = new TH1F("countsHist","Event Counts",10,0.0,10.0);
  countsHist->GetXaxis()->SetBinLabel(1,"total");
  countsHist->GetXaxis()->SetBinLabel(2,">=2#mu");
  countsHist->GetXaxis()->SetBinLabel(3,"VBFL");
  countsHist->GetXaxis()->SetBinLabel(4,"VBF");
  countsHist->GetXaxis()->SetBinLabel(5,"VBFT");
  countsHist->GetXaxis()->SetBinLabel(6,"Pt30");
  countsHist->GetXaxis()->SetBinLabel(7,"Pt50");
  countsHist->GetXaxis()->SetBinLabel(8,"Pt75");
  histMap.insert(make_pair("countsHist",countsHist));

  TH1F* cosThetaStarHist = new TH1F("cosThetaStar","cos(#theta^{*})",50,-1.,1.);
  histMap.insert(make_pair("cosThetaStar",cosThetaStarHist));
  TH1F* cosThetaStarCSHist = new TH1F("cosThetaStarCS","cos(#theta^{*}_{CS})",50,-1.,1.);
  histMap.insert(make_pair("cosThetaStarCS",cosThetaStarCSHist));

  TH1F* puJetIDSimpleDiscJet1Hist = new TH1F("puJetIDSimpleDiscJet1","PU Jet ID--Simple Discriminator Leading Jet",50,-1.,1.);
  histMap.insert(make_pair("puJetIDSimpleDiscJet1",puJetIDSimpleDiscJet1Hist));
  TH1F* puJetIDSimpleDiscJet2Hist = new TH1F("puJetIDSimpleDiscJet2","PU Jet ID--Simple Discriminator Sub-Leading Jet",50,-1.,1.);
  histMap.insert(make_pair("puJetIDSimpleDiscJet2",puJetIDSimpleDiscJet2Hist));
  TH1F* puJetIDSimpleDiscJet3Hist = new TH1F("puJetIDSimpleDiscJet3","PU Jet ID--Simple Discriminator 3rd Leading Jet",50,-1.,1.);
  histMap.insert(make_pair("puJetIDSimpleDiscJet3",puJetIDSimpleDiscJet3Hist));

  TH1F* puJetIDSimpleJet1Hist = new TH1F("puJetIDSimpleJet1","PU Jet ID--Simple Loose Leading Jet",2,0,2);
  histMap.insert(make_pair("puJetIDSimpleJet1",puJetIDSimpleJet1Hist));
  puJetIDSimpleJet1Hist->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet1Hist->GetXaxis()->SetBinLabel(2,"Pass");
  TH1F* puJetIDSimpleJet2Hist = new TH1F("puJetIDSimpleJet2","PU Jet ID--Simple Loose Sub-Leading Jet",2,-0,2);
  histMap.insert(make_pair("puJetIDSimpleJet2",puJetIDSimpleJet2Hist));
  puJetIDSimpleJet2Hist->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet2Hist->GetXaxis()->SetBinLabel(2,"Pass");
  TH1F* puJetIDSimpleJet3Hist = new TH1F("puJetIDSimpleJet3","PU Jet ID--Simple Loose 3rd Leading Jet",2,-0,2);
  histMap.insert(make_pair("puJetIDSimpleJet3",puJetIDSimpleJet3Hist));
  puJetIDSimpleJet3Hist->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet3Hist->GetXaxis()->SetBinLabel(2,"Pass");

  TH1F* BDTHistMuonOnly = new TH1F("BDTHistMuonOnly","BDT Discriminator",2000,-1,1);
  histMap.insert(make_pair("BDTHistMuonOnly",BDTHistMuonOnly));
  TH1F* likelihoodHistMuonOnly = new TH1F("likelihoodHistMuonOnly","Likelihood Discriminator",2000,-1,1);
  histMap.insert(make_pair("likelihoodHistMuonOnly",likelihoodHistMuonOnly));

  TH1F* BDTHistVBF = new TH1F("BDTHistVBF","BDT Discriminator",2000,-1,1);
  histMap.insert(make_pair("BDTHistVBF",BDTHistVBF));
  TH1F* likelihoodHistVBF = new TH1F("likelihoodHistVBF","Likelihood Discriminator",2000,-1,1);
  histMap.insert(make_pair("likelihoodHistVBF",likelihoodHistVBF));

  TH2F* BDTHistMuonOnlyVMass = new TH2F("BDTHistMuonOnlyVMass","BDT Discriminator",1600,0,400,2000,-1,1);
  histMap2D.insert(make_pair("BDTHistMuonOnlyVMass",BDTHistMuonOnlyVMass));
  TH2F* likelihoodHistMuonOnlyVMass = new TH2F("likelihoodHistMuonOnlyVMass","Likelihood Discriminator",1600,0,400,2000,-1,1);
  histMap2D.insert(make_pair("likelihoodHistMuonOnlyVMass",likelihoodHistMuonOnlyVMass));
  TH2F* BDTHistVBFVMass = new TH2F("BDTHistVBFVMass","BDT Discriminator",1600,0,400,2000,-1,1);
  histMap2D.insert(make_pair("BDTHistVBFVMass",BDTHistVBFVMass));
  TH2F* likelihoodHistVBFVMass = new TH2F("likelihoodHistVBFVMass","Likelihood Discriminator",1600,0,400,2000,-1,1);
  histMap2D.insert(make_pair("likelihoodHistVBFVMass",likelihoodHistVBFVMass));

  TH1F* relIsoMu1Hist = new TH1F("relIsoMu1","",1000,0,10.0);
  histMap.insert(make_pair("relIsoMu1",relIsoMu1Hist));
  TH1F* relIsoMu2Hist = new TH1F("relIsoMu2","",1000,0,10.0);
  histMap.insert(make_pair("relIsoMu2",relIsoMu2Hist));

  TH1F* nJetsHist = new TH1F("nJets","",11,0,11);
  histMap.insert(make_pair("nJets",nJetsHist));
  TH1F* htHist = new TH1F("ht","",200,0,2000);
  histMap.insert(make_pair("ht",htHist));
  TH1F* nJetsInRapidityGapHist = new TH1F("nJetsInRapidityGap","",11,0,11);
  histMap.insert(make_pair("nJetsInRapidityGap",nJetsInRapidityGapHist));
  TH1F* htInRapidityGapHist = new TH1F("htInRapidityGap","",200,0,2000);
  histMap.insert(make_pair("htInRapidityGap",htInRapidityGapHist));

  TH1F* nPUHist = new TH1F("nPU","",100,0,100);
  histMap.insert(make_pair("nPU",nPUHist));
  TH1F* nVtxHist = new TH1F("nVtx","",100,0,100);
  histMap.insert(make_pair("nVtx",nVtxHist));
  TH1F* metHist = new TH1F("met","",160,0,800);
  histMap.insert(make_pair("met",metHist));
  TH1F* weightHist = new TH1F("weight","",500,0,5.0);
  histMap.insert(make_pair("weight",weightHist));

  for(histMapIter = histMap.begin(); histMapIter != histMap.end(); histMapIter++)
  {
    histMapIter->second->Sumw2();
  }
  for(histMap2DIter = histMap2D.begin(); histMap2DIter != histMap2D.end(); histMap2DIter++)
  {
    histMap2DIter->second->Sumw2();
  }

  // Other Sets of Hists
  std::map<std::string,TH1F*> histMap4GeVWindow;
  std::map<std::string,TH1F*> histMapPtDiMu100;
  std::map<std::string,TH1F*> histMapVBFPresel;
  std::map<std::string,TH1F*> histMapIncPresel;
  std::map<std::string,TH1F*> histMapNotBlindWindow;
  std::map<std::string,TH1F*> histMapUpperControlRegion;
  std::map<std::string,TH1F*> histMapLowerControlRegion;
  for(histMapIter = histMap.begin(); histMapIter != histMap.end(); histMapIter++)
  {
    TH1F* tmp;
    tmp = (TH1F*) histMapIter->second->Clone();
    histMap4GeVWindow.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapPtDiMu100.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapVBFPresel.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapIncPresel.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapNotBlindWindow.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapUpperControlRegion.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapLowerControlRegion.insert(make_pair(histMapIter->first,tmp));
  }

  //////////////////////////
  //for MVA

  std::vector<std::string> mvaConfigNames;
  mvaConfigNames.push_back("inclusive.cfg");
  mvaConfigNames.push_back("vbf.cfg");
  MVA mva(mvaConfigNames,trainingTreeFileName);

  //////////////////////////
  //for PU reweighting

#ifdef PUREWEIGHT
  reweight::LumiReWeighting lumiWeights("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root","pileupDists/PileUpHist2012AB.root","pileup","pileup");
#endif

  const double SQRT2 = sqrt(2);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  unsigned nEvents = tree->GetEntries();
  cout << "nEvents: " << nEvents << endl;

  unsigned reportEach=1000;
  if (nEvents/1000>reportEach)
    reportEach = nEvents/1000;
  
  for(unsigned i=0; i<nEvents;i++)
  {
    if(i >= maxEvents)
      continue;
    tree->GetEvent(i);
    if (i % reportEach == 0) cout << "Event: " << i << endl;

    bool inBlindWindow = recoCandMass < maxBlind && recoCandMass > minBlind;
#ifdef BLIND
    if (inBlindWindow && isData)
        continue;
#endif

#ifdef PUREWEIGHT
    double weight = lumiWeights.weight(nPU);
#else
    double weight = 1.0;
#endif

    countsHist->Fill(0.0, weight);

    if (!isKinTight_2012(reco1) || !isKinTight_2012(reco2))
        continue;

    if (recoCandMass < minMmm || recoCandMass > maxMmm)
        continue;

    countsHist->Fill(1.0, weight);

    _MuonInfo muon1;
    _MuonInfo muon2;
    if(reco1.pt>reco2.pt)
    {
        muon1 = reco1;
        muon2 = reco2;
    }
    else
    {
        muon1 = reco2;
        muon2 = reco1;
    }
    mva.resetValues();

    mva.weight = weight;

    mva.ptMu1=muon1.pt;
    mva.ptMu2=muon2.pt;
    mva.etaMu1=muon1.eta;
    mva.etaMu2=muon2.eta;
    mva.deltaEtaMuons=fabs(muon1.eta-muon2.eta);
    mva.relIsoMu1 = getRelIso(muon1);
    mva.relIsoMu2 = getRelIso(muon2);

    mva.mDiMu = recoCandMass;
    mva.ptDiMu = recoCandPt;
    mva.yDiMu = recoCandY;

    //////////////////////////////////////////
    //Computing CosTheta*

    TLorentzVector pMuon1;
    TLorentzVector pMuon2;
    pMuon1.SetPtEtaPhiM(muon1.pt,muon1.eta,muon1.phi,0.105);
    pMuon2.SetPtEtaPhiM(muon2.pt,muon2.eta,muon2.phi,0.105);
    TLorentzVector diMuon = pMuon1+pMuon2;

    mva.deltaPhiMuons = pMuon1.DeltaPhi(pMuon2);
    mva.deltaRMuons = pMuon1.DeltaR(pMuon2);

    TLorentzVector starMuon1 = pMuon1;
    TLorentzVector starMuon2 = pMuon2;
    TVector3 boost = diMuon.BoostVector();
    starMuon1.Boost(-boost);
    starMuon2.Boost(-boost);

    //if (muon1.charge>0)
    //std::cout << "Run: " << eventInfo.run << " lumi: " << eventInfo.lumi << " event: " << eventInfo.event << std::endl;
    if ((int) (eventInfo.event) % 2 == 0)
    {
        TVector3 directionOfBoost = starMuon1.BoostVector();
        mva.cosThetaStar = directionOfBoost.Dot(diMuon.BoostVector()) / (directionOfBoost.Mag()*diMuon.BoostVector().Mag());
    }
    else
    {
        TVector3 directionOfBoost = starMuon2.BoostVector();
        mva.cosThetaStar = directionOfBoost.Dot(diMuon.BoostVector()) / (directionOfBoost.Mag()*diMuon.BoostVector().Mag());
    }

    //////////////////////////////////////////
    //Computing CosTheta* Collins-Soper

    //std::cout << "muon1 charge: " << muon1.charge << "muon2 charge: "<<muon2.charge << std::endl;
    if (muon1.charge != muon2.charge)
    {
      // p1 is lepton
      // p2 is anti-lepton
      float p1Plus=-1e15;
      float p2Plus=-1e15;
      float p1Minus=-1e15;
      float p2Minus=-1e15;
      if (muon1.charge < 0)
      {
        p1Plus  = (pMuon1.E()+pMuon1.Pz())/SQRT2;
        p1Minus = (pMuon1.E()-pMuon1.Pz())/SQRT2;
        p2Plus  = (pMuon2.E()+pMuon2.Pz())/SQRT2;
        p2Minus = (pMuon2.E()-pMuon2.Pz())/SQRT2;
      }
      else
      {
        p1Plus  = (pMuon2.E()+pMuon2.Pz())/SQRT2;
        p1Minus = (pMuon2.E()-pMuon2.Pz())/SQRT2;
        p2Plus  = (pMuon1.E()+pMuon1.Pz())/SQRT2;
        p2Minus = (pMuon1.E()-pMuon1.Pz())/SQRT2;
      }
      mva.cosThetaStarCS = diMuon.Pz()/fabs(diMuon.Pz()) * 
                               2*(p1Plus*p2Minus-p1Minus*p2Plus) / 
                               (diMuon.Mag()*sqrt(diMuon.Mag2()+diMuon.Pt()*diMuon.Pt()));
    }

    // Computing nVtx Valid
    //for(unsigned iVtx=0;iVtx<vertexInfo.nVertices;iVtx++)
    //{
    //  if(vertexInfo.isValid[i])
    //    mva.nVtx++;
    //}

    //////////////////////////////////////////
    // Filling Hists

    mDiMu->Fill(mva.mDiMu, weight);
    yDiMu->Fill(mva.yDiMu, weight);
    ptDiMu->Fill(mva.ptDiMu, weight);
    yVptDiMu->Fill(mva.ptDiMu,fabs(mva.yDiMu), weight);
    yVmDiMu->Fill(mva.mDiMu,fabs(mva.yDiMu), weight);
    ptVmDiMu->Fill(mva.mDiMu,mva.ptDiMu, weight);
    phiVmDiMu->Fill(mva.mDiMu,recoCandPhi, weight);
    ptMu1->Fill(mva.ptMu1, weight);
    ptMu2->Fill(mva.ptMu2, weight);
    etaMu1->Fill(mva.etaMu1, weight);
    etaMu2->Fill(mva.etaMu2, weight);
    cosThetaStarHist->Fill(mva.cosThetaStar, weight);
    cosThetaStarCSHist->Fill(mva.cosThetaStarCS, weight);

    deltaPhiMuonsHist->Fill(mva.deltaPhiMuons, weight);
    deltaEtaMuonsHist->Fill(mva.deltaEtaMuons, weight);
    deltaRMuonsHist->Fill(mva.deltaRMuons, weight);

    relIsoMu1Hist->Fill(mva.relIsoMu1, weight);
    relIsoMu2Hist->Fill(mva.relIsoMu2, weight);

    nPUHist->Fill(nPU, weight);
    nVtxHist->Fill(mva.nVtx, weight);
    metHist->Fill(met.pt, weight);
    weightHist->Fill(weight);

    // Jet Part
    for(unsigned iJet=0; (iJet < jets.nPFjets && iJet < 10);iJet++)
    {
        if(jets.pfJetPt[iJet] > 30.0)
        {
          mva.nJets++;
          mva.ht += jets.pfJetPt[iJet];
        }
    }

    bool goodJets = false;
    if(jets.nPFjets>=2 && jets.pfJetPt[0]>30.0 && jets.pfJetPt[1]>30.0)
        goodJets = true;

    if(goodJets)
    {
      TLorentzVector pJet1;
      TLorentzVector pJet2;
      pJet1.SetXYZM(jets.pfJetPx[0],jets.pfJetPy[0],jets.pfJetPz[0],jets.pfJetM[0]);
      pJet2.SetXYZM(jets.pfJetPx[1],jets.pfJetPy[1],jets.pfJetPz[1],jets.pfJetM[1]);
      TLorentzVector diJet = pJet1+pJet2;

      double dEtaJets = fabs(jets.pfJetEta[0]-jets.pfJetEta[1]);
      double etaJetProduct = jets.pfJetEta[0]*jets.pfJetEta[1];
      mva.deltaPhiJets = pJet1.DeltaPhi(pJet2);
      mva.deltaRJets = pJet1.DeltaR(pJet2);

      // Seeing if there are jets in the rapidity gap
      float etaMax = jets.pfJetEta[0];
      float etaMin = 9999999.0;
      if(etaMax < jets.pfJetEta[1])
      {
          etaMax = jets.pfJetEta[1];
          etaMin = jets.pfJetEta[0];
      }
      else
      {
          etaMin = jets.pfJetEta[1];
      }
      bool jetInRapidityGap=false;
      for(unsigned iJet=2; (iJet < jets.nPFjets && iJet < 10);iJet++)
      {
        if(jets.pfJetPt[iJet] > 30.0)
        {
          if(jets.pfJetEta[iJet] < etaMax && jets.pfJetEta[iJet] > etaMin)
          {
            jetInRapidityGap = true;
            mva.nJetsInRapidityGap++;
            mva.htInRapidityGap += jets.pfJetPt[iJet];
          }
        }
      }

#ifdef JETPUID
      for(unsigned iJet=0; (iJet < jets.nPFjets && iJet < 10);iJet++)
      {
        if(jets.pfJetPt[iJet]>30.0)
        {
          if (iJet==0)
            puJetIDSimpleDiscJet1Hist->Fill(puJetSimpleDisc[iJet], weight);
          else if (iJet==1)
            puJetIDSimpleDiscJet2Hist->Fill(puJetSimpleDisc[iJet], weight);
          else if (iJet==2)
            puJetIDSimpleDiscJet3Hist->Fill(puJetSimpleDisc[iJet], weight);

          if (iJet==0)
            puJetIDSimpleJet1Hist->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose), weight);
          else if (iJet==1)
            puJetIDSimpleJet2Hist->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose), weight);
          else if (iJet==2)
            puJetIDSimpleJet3Hist->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose), weight);
        }
      }
#endif

      mva.mDiJet = diJet.M();
      mva.yDiJet = diJet.Rapidity();
      mva.ptDiJet = diJet.Pt();
      mva.ptJet1 = pJet1.Pt();
      mva.ptJet2 = pJet2.Pt();
      mva.etaJet1 = pJet1.Eta();
      mva.etaJet2 = pJet2.Eta();
      mva.productEtaJets = etaJetProduct;
      mva.deltaEtaJets = dEtaJets;

      mDiJet->Fill(mva.mDiJet, weight);
      ptJet1->Fill(mva.ptJet1, weight);
      ptJet2->Fill(mva.ptJet2, weight);
      etaJet1->Fill(mva.etaJet1, weight);
      etaJet2->Fill(mva.etaJet2, weight);
      deltaEtaJetsHist->Fill(mva.deltaEtaJets, weight);
      deltaPhiJetsHist->Fill(mva.deltaPhiJets, weight);
      deltaRJetsHist->Fill(mva.deltaRJets, weight);
      nJetsInRapidityGapHist->Fill(mva.nJetsInRapidityGap, weight);
      htInRapidityGapHist->Fill(mva.htInRapidityGap, weight);
      nJetsHist->Fill(mva.nJets, weight);
      htHist->Fill(mva.ht, weight);
    }
  
//HIG-12-007 PAS H->tautau
//The VBF category requires at least two jets with pT > 30 GeV/c, |η1 − η2 | > 4.0 and
//η1 · η2 < 0 (with η1 the pseudorapidty of the leading jet and η2 the pseudorapidity
//of the subleading jet), and a di-jet invariant mass m12 > 400 GeV/c2 , with no other
//jet with pT > 30 GeV/c in the rapidity region between the two jets.

    mva.writeEvent();

    bool vbfPreselection = mva.mDiJet>300.0 && mva.deltaEtaJets>3.0 && mva.productEtaJets<0.0 && mva.nJetsInRapidityGap == 0;
    //if(vbfPreselection)
    //  std::cout << "VBF Preselected!!";

    if(!vbfPreselection)
    {
      BDTHistMuonOnly->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
      likelihoodHistMuonOnly->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      BDTHistMuonOnlyVMass->Fill(recoCandMass, mva.getMVA("inclusive.cfg","BDT"), weight);
      likelihoodHistMuonOnlyVMass->Fill(recoCandMass, mva.getMVA("inclusive.cfg","Likelihood"), weight);
    }
    else
    {
      BDTHistVBFVMass->Fill(recoCandMass, mva.getMVA("vbf.cfg","BDT"), weight);
      likelihoodHistVBFVMass->Fill(recoCandMass, mva.getMVA("vbf.cfg","Likelihood"), weight);
      BDTHistVBF->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
      likelihoodHistVBF->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
    }

    //4 GeV Window Plots
    if (mva.mDiMu < 127.0 && mva.mDiMu > 123.0)
    {
      histMap4GeVWindow["mDiMu"]->Fill(mva.mDiMu, weight);
      histMap4GeVWindow["yDiMu"]->Fill(mva.yDiMu, weight);
      histMap4GeVWindow["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMap4GeVWindow["ptMu1"]->Fill(mva.ptMu1, weight);
      histMap4GeVWindow["ptMu2"]->Fill(mva.ptMu2, weight);
      histMap4GeVWindow["etaMu1"]->Fill(mva.etaMu1, weight);
      histMap4GeVWindow["etaMu2"]->Fill(mva.etaMu2, weight);
      histMap4GeVWindow["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMap4GeVWindow["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMap4GeVWindow["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMap4GeVWindow["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMap4GeVWindow["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMap4GeVWindow["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMap4GeVWindow["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMap4GeVWindow["nPU"]->Fill(nPU, weight);
      histMap4GeVWindow["nVtx"]->Fill(mva.nVtx, weight);
      histMap4GeVWindow["met"]->Fill(met.pt, weight);

      histMap4GeVWindow["mDiJet"]->Fill(mva.mDiJet, weight);
      histMap4GeVWindow["ptJet1"]->Fill(mva.ptJet1, weight);
      histMap4GeVWindow["ptJet2"]->Fill(mva.ptJet2, weight);
      histMap4GeVWindow["etaJet1"]->Fill(mva.etaJet1, weight);
      histMap4GeVWindow["etaJet2"]->Fill(mva.etaJet2, weight);
      histMap4GeVWindow["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMap4GeVWindow["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMap4GeVWindow["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMap4GeVWindow["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMap4GeVWindow["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMap4GeVWindow["nJets"]->Fill(mva.nJets, weight);
      histMap4GeVWindow["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMap4GeVWindow["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMap4GeVWindow["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMap4GeVWindow["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMap4GeVWindow["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }

    //DiMu Pt > 100 Plots
    if (mva.ptDiMu > 100.0)
    {
      histMapPtDiMu100["mDiMu"]->Fill(mva.mDiMu, weight);
      histMapPtDiMu100["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapPtDiMu100["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapPtDiMu100["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapPtDiMu100["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapPtDiMu100["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapPtDiMu100["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapPtDiMu100["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapPtDiMu100["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapPtDiMu100["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapPtDiMu100["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapPtDiMu100["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapPtDiMu100["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapPtDiMu100["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapPtDiMu100["nPU"]->Fill(nPU, weight);
      histMapPtDiMu100["nVtx"]->Fill(mva.nVtx, weight);
      histMapPtDiMu100["met"]->Fill(met.pt, weight);

      histMapPtDiMu100["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapPtDiMu100["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapPtDiMu100["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapPtDiMu100["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapPtDiMu100["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapPtDiMu100["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapPtDiMu100["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapPtDiMu100["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapPtDiMu100["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapPtDiMu100["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapPtDiMu100["nJets"]->Fill(mva.nJets, weight);
      histMapPtDiMu100["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMapPtDiMu100["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapPtDiMu100["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapPtDiMu100["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapPtDiMu100["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }

    //VBF Preselected Plots
    if (vbfPreselection)
    {
      histMapVBFPresel["mDiMu"]->Fill(mva.mDiMu, weight);
      histMapVBFPresel["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapVBFPresel["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapVBFPresel["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapVBFPresel["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapVBFPresel["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapVBFPresel["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapVBFPresel["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapVBFPresel["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapVBFPresel["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapVBFPresel["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapVBFPresel["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapVBFPresel["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapVBFPresel["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapVBFPresel["nPU"]->Fill(nPU, weight);
      histMapVBFPresel["nVtx"]->Fill(mva.nVtx, weight);
      histMapVBFPresel["met"]->Fill(met.pt, weight);

      histMapVBFPresel["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapVBFPresel["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapVBFPresel["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapVBFPresel["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapVBFPresel["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapVBFPresel["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapVBFPresel["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapVBFPresel["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapVBFPresel["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapVBFPresel["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapVBFPresel["nJets"]->Fill(mva.nJets, weight);
      histMapVBFPresel["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMapVBFPresel["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapVBFPresel["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapVBFPresel["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapVBFPresel["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }

    //Not VBF Preselected Plots
    if (!vbfPreselection)
    {
      histMapIncPresel["mDiMu"]->Fill(mva.mDiMu, weight);
      histMapIncPresel["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapIncPresel["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapIncPresel["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapIncPresel["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapIncPresel["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapIncPresel["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapIncPresel["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapIncPresel["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapIncPresel["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapIncPresel["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapIncPresel["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapIncPresel["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapIncPresel["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapIncPresel["nPU"]->Fill(nPU, weight);
      histMapIncPresel["nVtx"]->Fill(mva.nVtx, weight);
      histMapIncPresel["met"]->Fill(met.pt, weight);

      histMapIncPresel["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapIncPresel["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapIncPresel["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapIncPresel["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapIncPresel["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapIncPresel["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapIncPresel["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapIncPresel["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapIncPresel["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapIncPresel["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapIncPresel["nJets"]->Fill(mva.nJets, weight);
      histMapIncPresel["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMapIncPresel["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapIncPresel["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapIncPresel["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapIncPresel["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }

    //Not In Blinding Window Plots
    if (!inBlindWindow)
    {
      histMapNotBlindWindow["mDiMu"]->Fill(mva.mDiMu, weight);
      histMapNotBlindWindow["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapNotBlindWindow["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapNotBlindWindow["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapNotBlindWindow["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapNotBlindWindow["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapNotBlindWindow["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapNotBlindWindow["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapNotBlindWindow["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapNotBlindWindow["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapNotBlindWindow["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapNotBlindWindow["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapNotBlindWindow["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapNotBlindWindow["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapNotBlindWindow["nPU"]->Fill(nPU, weight);
      histMapNotBlindWindow["nVtx"]->Fill(mva.nVtx, weight);
      histMapNotBlindWindow["met"]->Fill(met.pt, weight);

      histMapNotBlindWindow["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapNotBlindWindow["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapNotBlindWindow["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapNotBlindWindow["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapNotBlindWindow["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapNotBlindWindow["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapNotBlindWindow["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapNotBlindWindow["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapNotBlindWindow["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapNotBlindWindow["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapNotBlindWindow["nJets"]->Fill(mva.nJets, weight);
      histMapNotBlindWindow["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMapNotBlindWindow["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapNotBlindWindow["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapNotBlindWindow["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapNotBlindWindow["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }

    //Upper Control Region Plots
    if (recoCandMass>maxBlind)
    {
      histMapUpperControlRegion["mDiMu"]->Fill(mva.mDiMu, weight);
      histMapUpperControlRegion["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapUpperControlRegion["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapUpperControlRegion["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapUpperControlRegion["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapUpperControlRegion["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapUpperControlRegion["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapUpperControlRegion["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapUpperControlRegion["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapUpperControlRegion["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapUpperControlRegion["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapUpperControlRegion["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapUpperControlRegion["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapUpperControlRegion["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapUpperControlRegion["nPU"]->Fill(nPU, weight);
      histMapUpperControlRegion["nVtx"]->Fill(mva.nVtx, weight);
      histMapUpperControlRegion["met"]->Fill(met.pt, weight);

      histMapUpperControlRegion["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapUpperControlRegion["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapUpperControlRegion["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapUpperControlRegion["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapUpperControlRegion["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapUpperControlRegion["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapUpperControlRegion["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapUpperControlRegion["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapUpperControlRegion["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapUpperControlRegion["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapUpperControlRegion["nJets"]->Fill(mva.nJets, weight);
      histMapUpperControlRegion["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMapUpperControlRegion["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapUpperControlRegion["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapUpperControlRegion["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapUpperControlRegion["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }
    //Lower Control Region Plots
    if (recoCandMass<minBlind)
    {
      histMapLowerControlRegion["mDiMu"]->Fill(mva.mDiMu, weight);
      histMapLowerControlRegion["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapLowerControlRegion["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapLowerControlRegion["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapLowerControlRegion["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapLowerControlRegion["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapLowerControlRegion["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapLowerControlRegion["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapLowerControlRegion["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapLowerControlRegion["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapLowerControlRegion["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapLowerControlRegion["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapLowerControlRegion["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapLowerControlRegion["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapLowerControlRegion["nPU"]->Fill(nPU, weight);
      histMapLowerControlRegion["nVtx"]->Fill(mva.nVtx, weight);
      histMapLowerControlRegion["met"]->Fill(met.pt, weight);

      histMapLowerControlRegion["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapLowerControlRegion["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapLowerControlRegion["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapLowerControlRegion["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapLowerControlRegion["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapLowerControlRegion["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapLowerControlRegion["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapLowerControlRegion["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapLowerControlRegion["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapLowerControlRegion["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapLowerControlRegion["nJets"]->Fill(mva.nJets, weight);
      histMapLowerControlRegion["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMapLowerControlRegion["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapLowerControlRegion["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapLowerControlRegion["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapLowerControlRegion["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }

  }// end event loop

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  outFile->cd();

  for(histMapIter = histMap.begin(); histMapIter != histMap.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  for(histMap2DIter = histMap2D.begin(); histMap2DIter != histMap2D.end(); histMap2DIter++)
  {
    histMap2DIter->second->Write();
  }

  TDirectory* dir4GeVWindow = outFile->mkdir("4GeVWindow");
  dir4GeVWindow->cd();
  for(histMapIter = histMap4GeVWindow.begin(); histMapIter != histMap4GeVWindow.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirPtDiMu100 = outFile->mkdir("PtDiMu100");
  dirPtDiMu100->cd();
  for(histMapIter = histMapPtDiMu100.begin(); histMapIter != histMapPtDiMu100.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirVBFPresel = outFile->mkdir("VBFPresel");
  dirVBFPresel->cd();
  for(histMapIter = histMapVBFPresel.begin(); histMapIter != histMapVBFPresel.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirIncPresel = outFile->mkdir("IncPresel");
  dirIncPresel->cd();
  for(histMapIter = histMapIncPresel.begin(); histMapIter != histMapIncPresel.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirNotBlindWindow = outFile->mkdir("NotBlindWindow");
  dirNotBlindWindow->cd();
  for(histMapIter = histMapNotBlindWindow.begin(); histMapIter != histMapNotBlindWindow.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirUpperControlRegion = outFile->mkdir("UpperControlRegion");
  dirUpperControlRegion->cd();
  for(histMapIter = histMapUpperControlRegion.begin(); histMapIter != histMapUpperControlRegion.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirLowerControlRegion = outFile->mkdir("LowerControlRegion");
  dirLowerControlRegion->cd();
  for(histMapIter = histMapLowerControlRegion.begin(); histMapIter != histMapLowerControlRegion.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  cout << "analyzer done." << endl << endl;
  return 0;
}
