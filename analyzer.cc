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

#include "boost/program_options.hpp"
#include "boost/regex.hpp"

#include <limits.h>

#define JETPUID

using namespace std;
using namespace boost;

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  const char* optionIntro = "H->MuMu Analyzer\n\nUsage: ./analyzer [--help] [--train] [--maxEvents N] <inputFileName.root> <outputFileName.root>\n\nAllowed Options";
  program_options::options_description optionDesc(optionIntro);
  optionDesc.add_options()
      ("help,h", "Produce Help Message")
      ("trainingTree,t",program_options::value<string>(), "Create Training Tree File with filename")
      ("filenames",program_options::value<vector<string> >(), "Input & Output File Names")
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
  
  std::string inputFileName;
  std::string outputFileName;
  if (optionMap.count("filenames")>0)
  {
     vector<string> filenames = optionMap["filenames"].as<vector<string> >();
     if(filenames.size()>2)
     {
       cout << "Error: Extra filenames on command line, exiting." << endl;
       return 1;
     }
     if(filenames.size()<2)
     {
       cout << "Error: Need both input file and output file names, exiting." << endl;
       return 1;
     }
     inputFileName = filenames[0];
     outputFileName = filenames[1];
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

  cout << "Input File Name: " << inputFileName << endl;
  cout << "Output File Name: " << outputFileName << endl;

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
    //bool tmpIsData = regex_match(inputFileName,re);
    //if(tmpIsData)
    if(regex_search(inputFileName,re))
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
  tree->AddFile(inputFileName.c_str());

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

  float recoCandMass, recoCandPt, recoCandY;

  tree->SetBranchAddress("recoCandMass", &recoCandMass);
  tree->SetBranchAddress("recoCandPt"  , &recoCandPt );
  tree->SetBranchAddress("recoCandY"  , &recoCandY );

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

  // Small Mass Window Hists
  std::map<std::string,TH1F*> histMap4GeVWindow;
  std::map<std::string,TH1F*> histMapPtDiMu100;
  std::map<std::string,TH1F*> histMapVBFPresel;
  std::map<std::string,TH1F*> histMapIncPresel;
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
  }

  //////////////////////////
  //for MVA

  std::vector<std::string> mvaConfigNames;
  mvaConfigNames.push_back("inclusive.cfg");
  mvaConfigNames.push_back("vbf.cfg");
  MVA mva(mvaConfigNames,trainingTreeFileName);


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

#ifdef BLIND
    bool inBlindWindow = recoCandMass < 145.0 && recoCandMass > 105.0;
    if (inBlindWindow && isData)
        continue;
#endif

    countsHist->Fill(0.0);

    if (!isKinTight_2012(reco1) || !isKinTight_2012(reco2))
        continue;

    if (recoCandMass < minMmm || recoCandMass > maxMmm)
        continue;

    countsHist->Fill(1.0);

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
    if ((int) (1000 * muon1.pt) % 2 == 0)
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
    // Filling Hists

    mDiMu->Fill(mva.mDiMu);
    yDiMu->Fill(mva.yDiMu);
    ptDiMu->Fill(mva.ptDiMu);
    yVptDiMu->Fill(mva.ptDiMu,fabs(mva.yDiMu));
    ptMu1->Fill(mva.ptMu1);
    ptMu2->Fill(mva.ptMu2);
    etaMu1->Fill(mva.etaMu1);
    etaMu2->Fill(mva.etaMu2);
    cosThetaStarHist->Fill(mva.cosThetaStar);

    deltaPhiMuonsHist->Fill(mva.deltaPhiMuons);
    deltaEtaMuonsHist->Fill(mva.deltaEtaMuons);
    deltaRMuonsHist->Fill(mva.deltaRMuons);

    relIsoMu1Hist->Fill(mva.relIsoMu1);
    relIsoMu2Hist->Fill(mva.relIsoMu2);

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
            puJetIDSimpleDiscJet1Hist->Fill(puJetSimpleDisc[iJet]);
          else if (iJet==1)
            puJetIDSimpleDiscJet2Hist->Fill(puJetSimpleDisc[iJet]);
          else if (iJet==2)
            puJetIDSimpleDiscJet3Hist->Fill(puJetSimpleDisc[iJet]);

          if (iJet==0)
            puJetIDSimpleJet1Hist->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose));
          else if (iJet==1)
            puJetIDSimpleJet2Hist->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose));
          else if (iJet==2)
            puJetIDSimpleJet3Hist->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose));
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

      mDiJet->Fill(mva.mDiJet);
      ptJet1->Fill(mva.ptJet1);
      ptJet2->Fill(mva.ptJet2);
      etaJet1->Fill(mva.etaJet1);
      etaJet2->Fill(mva.etaJet1);
      deltaEtaJetsHist->Fill(mva.deltaEtaJets);
      deltaPhiJetsHist->Fill(mva.deltaPhiJets);
      deltaRJetsHist->Fill(mva.deltaRJets);
      nJetsInRapidityGapHist->Fill(mva.nJetsInRapidityGap);
      htInRapidityGapHist->Fill(mva.htInRapidityGap);
      nJetsHist->Fill(mva.nJets);
      htHist->Fill(mva.ht);
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
      BDTHistMuonOnly->Fill(mva.getMVA("inclusive.cfg","BDT"));
      likelihoodHistMuonOnly->Fill(mva.getMVA("inclusive.cfg","Likelihood"));
      BDTHistMuonOnlyVMass->Fill(recoCandMass, mva.getMVA("inclusive.cfg","BDT"));
      likelihoodHistMuonOnlyVMass->Fill(recoCandMass, mva.getMVA("inclusive.cfg","Likelihood"));
    }
    else
    {
      BDTHistVBFVMass->Fill(recoCandMass, mva.getMVA("vbf.cfg","BDT"));
      likelihoodHistVBFVMass->Fill(recoCandMass, mva.getMVA("vbf.cfg","Likelihood"));
      BDTHistVBF->Fill(mva.getMVA("vbf.cfg","BDT"));
      likelihoodHistVBF->Fill(mva.getMVA("vbf.cfg","Likelihood"));
    }

    //4 GeV Window Plots
    if (mva.mDiMu < 127.0 && mva.mDiMu > 123.0)
    {
      histMap4GeVWindow["mDiMu"]->Fill(mva.mDiMu);
      histMap4GeVWindow["yDiMu"]->Fill(mva.yDiMu);
      histMap4GeVWindow["ptDiMu"]->Fill(mva.ptDiMu);
      histMap4GeVWindow["ptMu1"]->Fill(mva.ptMu1);
      histMap4GeVWindow["ptMu2"]->Fill(mva.ptMu2);
      histMap4GeVWindow["etaMu1"]->Fill(mva.etaMu1);
      histMap4GeVWindow["etaMu2"]->Fill(mva.etaMu2);
      histMap4GeVWindow["cosThetaStar"]->Fill(mva.cosThetaStar);
      histMap4GeVWindow["deltaPhiMuons"]->Fill(mva.deltaPhiMuons);
      histMap4GeVWindow["deltaEtaMuons"]->Fill(mva.deltaEtaMuons);
      histMap4GeVWindow["deltaRMuons"]->Fill(mva.deltaRMuons);
      histMap4GeVWindow["relIsoMu1"]->Fill(mva.relIsoMu1);
      histMap4GeVWindow["relIsoMu2"]->Fill(mva.relIsoMu2);

      histMap4GeVWindow["mDiJet"]->Fill(mva.mDiJet);
      histMap4GeVWindow["ptJet1"]->Fill(mva.ptJet1);
      histMap4GeVWindow["ptJet2"]->Fill(mva.ptJet2);
      histMap4GeVWindow["etaJet1"]->Fill(mva.etaJet1);
      histMap4GeVWindow["etaJet2"]->Fill(mva.etaJet1);
      histMap4GeVWindow["deltaEtaJets"]->Fill(mva.deltaEtaJets);
      histMap4GeVWindow["deltaPhiJets"]->Fill(mva.deltaPhiJets);
      histMap4GeVWindow["deltaRJets"]->Fill(mva.deltaRJets);
      histMap4GeVWindow["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap);
      histMap4GeVWindow["htInRapidityGap"]->Fill(mva.htInRapidityGap);
      histMap4GeVWindow["nJets"]->Fill(mva.nJets);
      histMap4GeVWindow["ht"]->Fill(mva.ht);
    }

    //DiMu Pt > 100 Plots
    if (mva.ptDiMu > 100.0)
    {
      histMapPtDiMu100["mDiMu"]->Fill(mva.mDiMu);
      histMapPtDiMu100["yDiMu"]->Fill(mva.yDiMu);
      histMapPtDiMu100["ptDiMu"]->Fill(mva.ptDiMu);
      histMapPtDiMu100["ptMu1"]->Fill(mva.ptMu1);
      histMapPtDiMu100["ptMu2"]->Fill(mva.ptMu2);
      histMapPtDiMu100["etaMu1"]->Fill(mva.etaMu1);
      histMapPtDiMu100["etaMu2"]->Fill(mva.etaMu2);
      histMapPtDiMu100["cosThetaStar"]->Fill(mva.cosThetaStar);
      histMapPtDiMu100["deltaPhiMuons"]->Fill(mva.deltaPhiMuons);
      histMapPtDiMu100["deltaEtaMuons"]->Fill(mva.deltaEtaMuons);
      histMapPtDiMu100["deltaRMuons"]->Fill(mva.deltaRMuons);
      histMapPtDiMu100["relIsoMu1"]->Fill(mva.relIsoMu1);
      histMapPtDiMu100["relIsoMu2"]->Fill(mva.relIsoMu2);

      histMapPtDiMu100["mDiJet"]->Fill(mva.mDiJet);
      histMapPtDiMu100["ptJet1"]->Fill(mva.ptJet1);
      histMapPtDiMu100["ptJet2"]->Fill(mva.ptJet2);
      histMapPtDiMu100["etaJet1"]->Fill(mva.etaJet1);
      histMapPtDiMu100["etaJet2"]->Fill(mva.etaJet1);
      histMapPtDiMu100["deltaEtaJets"]->Fill(mva.deltaEtaJets);
      histMapPtDiMu100["deltaPhiJets"]->Fill(mva.deltaPhiJets);
      histMapPtDiMu100["deltaRJets"]->Fill(mva.deltaRJets);
      histMapPtDiMu100["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap);
      histMapPtDiMu100["htInRapidityGap"]->Fill(mva.htInRapidityGap);
      histMapPtDiMu100["nJets"]->Fill(mva.nJets);
      histMapPtDiMu100["ht"]->Fill(mva.ht);
    }

    //VBF Preselected Plots
    if (vbfPreselection)
    {
      histMapVBFPresel["mDiMu"]->Fill(mva.mDiMu);
      histMapVBFPresel["yDiMu"]->Fill(mva.yDiMu);
      histMapVBFPresel["ptDiMu"]->Fill(mva.ptDiMu);
      histMapVBFPresel["ptMu1"]->Fill(mva.ptMu1);
      histMapVBFPresel["ptMu2"]->Fill(mva.ptMu2);
      histMapVBFPresel["etaMu1"]->Fill(mva.etaMu1);
      histMapVBFPresel["etaMu2"]->Fill(mva.etaMu2);
      histMapVBFPresel["cosThetaStar"]->Fill(mva.cosThetaStar);
      histMapVBFPresel["deltaPhiMuons"]->Fill(mva.deltaPhiMuons);
      histMapVBFPresel["deltaEtaMuons"]->Fill(mva.deltaEtaMuons);
      histMapVBFPresel["deltaRMuons"]->Fill(mva.deltaRMuons);
      histMapVBFPresel["relIsoMu1"]->Fill(mva.relIsoMu1);
      histMapVBFPresel["relIsoMu2"]->Fill(mva.relIsoMu2);

      histMapVBFPresel["mDiJet"]->Fill(mva.mDiJet);
      histMapVBFPresel["ptJet1"]->Fill(mva.ptJet1);
      histMapVBFPresel["ptJet2"]->Fill(mva.ptJet2);
      histMapVBFPresel["etaJet1"]->Fill(mva.etaJet1);
      histMapVBFPresel["etaJet2"]->Fill(mva.etaJet1);
      histMapVBFPresel["deltaEtaJets"]->Fill(mva.deltaEtaJets);
      histMapVBFPresel["deltaPhiJets"]->Fill(mva.deltaPhiJets);
      histMapVBFPresel["deltaRJets"]->Fill(mva.deltaRJets);
      histMapVBFPresel["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap);
      histMapVBFPresel["htInRapidityGap"]->Fill(mva.htInRapidityGap);
      histMapVBFPresel["nJets"]->Fill(mva.nJets);
      histMapVBFPresel["ht"]->Fill(mva.ht);
    }

    //Not VBF Preselected Plots
    if (!vbfPreselection)
    {
      histMapIncPresel["mDiMu"]->Fill(mva.mDiMu);
      histMapIncPresel["yDiMu"]->Fill(mva.yDiMu);
      histMapIncPresel["ptDiMu"]->Fill(mva.ptDiMu);
      histMapIncPresel["ptMu1"]->Fill(mva.ptMu1);
      histMapIncPresel["ptMu2"]->Fill(mva.ptMu2);
      histMapIncPresel["etaMu1"]->Fill(mva.etaMu1);
      histMapIncPresel["etaMu2"]->Fill(mva.etaMu2);
      histMapIncPresel["cosThetaStar"]->Fill(mva.cosThetaStar);
      histMapIncPresel["deltaPhiMuons"]->Fill(mva.deltaPhiMuons);
      histMapIncPresel["deltaEtaMuons"]->Fill(mva.deltaEtaMuons);
      histMapIncPresel["deltaRMuons"]->Fill(mva.deltaRMuons);
      histMapIncPresel["relIsoMu1"]->Fill(mva.relIsoMu1);
      histMapIncPresel["relIsoMu2"]->Fill(mva.relIsoMu2);

      histMapIncPresel["mDiJet"]->Fill(mva.mDiJet);
      histMapIncPresel["ptJet1"]->Fill(mva.ptJet1);
      histMapIncPresel["ptJet2"]->Fill(mva.ptJet2);
      histMapIncPresel["etaJet1"]->Fill(mva.etaJet1);
      histMapIncPresel["etaJet2"]->Fill(mva.etaJet1);
      histMapIncPresel["deltaEtaJets"]->Fill(mva.deltaEtaJets);
      histMapIncPresel["deltaPhiJets"]->Fill(mva.deltaPhiJets);
      histMapIncPresel["deltaRJets"]->Fill(mva.deltaRJets);
      histMapIncPresel["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap);
      histMapIncPresel["htInRapidityGap"]->Fill(mva.htInRapidityGap);
      histMapIncPresel["nJets"]->Fill(mva.nJets);
      histMapIncPresel["ht"]->Fill(mva.ht);
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

  cout << "analyzer done." << endl << endl;
  return 0;
}
