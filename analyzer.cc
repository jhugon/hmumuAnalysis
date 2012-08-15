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
#include "Helpers.h"

//#define JETPUID

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 3) 
  {
	std::cout << "testReader requires >=2 argument:" << std::endl
	<< "analyzer <outfilename.root> <infilenames.root...>" << std::endl;
	return 1;
  }

  float minMmm = 100.0;
  float maxMmm = 200.0;
  //minMmm = 123.0;
  //maxMmm = 127.0;

  TChain * tree = new TChain("tree");

  TFile * outFile = new TFile(argv[1],"RECREATE");

  long maxEvents = 100000000;
  //maxEvents = 20000;

  for(unsigned iarg=2; iarg<argc;iarg++)
  {
    tree->AddFile(argv[iarg]);
    cout << "Running on file " << argv[iarg] << endl;
  }
  cout << "Output file " << argv[1] << endl;

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

  TH1F* mDiMu = new TH1F("mDiMu","DiMuon Mass",1600,0,400);
  TH1F* mDiMuVBFM = new TH1F("mDiMuVBFM","DiMuon Mass after VBF Selection",1600,0,400);
  TH1F* mDiMuVBFL = new TH1F("mDiMuVBFL","DiMuon Mass after VBFLoose Selection",1600,0,400);
  TH1F* mDiMuVBFT = new TH1F("mDiMuVBFT","DiMuon Mass after VBFTight Selection",1600,0,400);
  TH1F* mDiMuVBFVL = new TH1F("mDiMuVBFVL","DiMuon Mass after VBFVeryLoose Selection",1600,0,400);
  TH1F* mDiMuPtL30 = new TH1F("mDiMuPtL30","DiMuon Mass after p_T^{#mu#mu}>30 GeV Selection",1600,0,400);
  TH1F* mDiMuPt30to50 = new TH1F("mDiMuPt30to50","DiMuon Mass after p_T^{#mu#mu}>30 GeV Selection",1600,0,400);
  TH1F* mDiMuPt50to75 = new TH1F("mDiMuPt50to75","DiMuon Mass after p_T^{#mu#mu}>50 GeV Selection",1600,0,400);
  TH1F* mDiMuPt75to125 = new TH1F("mDiMuPt75to125","DiMuon Mass after p_T^{#mu#mu}>50 GeV Selection",1600,0,400);
  TH1F* mDiMuPt125 = new TH1F("mDiMuPt125","DiMuon Mass after p_T^{#mu#mu}>75 GeV Selection",1600,0,400);

  TH1F* mDiMuEta11 = new TH1F("mDiMuEta11","DiMuon Mass",1600,0,400);
  TH1F* mDiMuEta12 = new TH1F("mDiMuEta12","DiMuon Mass",1600,0,400);
  TH1F* mDiMuEta22 = new TH1F("mDiMuEta22","DiMuon Mass",1600,0,400);

  TH1F* mDiJet = new TH1F("mDiJet","DiJet Mass",500,0,2000);

  TH1F* ptDiMu = new TH1F("ptDiMu","DiMuon Pt",250,0,500);
  TH1F* ptDiMuVBFM = new TH1F("ptDiMuVBFM","DiMuon Pt after VBF Selection",250,0,500);
  TH1F* ptDiMuVBFL = new TH1F("ptDiMuVBFL","DiMuon Pt after VBFLoose Selection",250,0,500);

  TH1F* yDiMu = new TH1F("yDiMu","DiMuon Rapidity",100,-4,4);
  TH1F* yDiMuVBFM = new TH1F("yDiMuVBFM","DiMuon Rapidity after VBF Selection",100,-4,4);
  TH1F* yDiMuVBFL = new TH1F("yDiMuVBFL","DiMuon Rapidity after VBFLoose Selection",100,-4,4);
  TH1F* yDiMuPt30 = new TH1F("yDiMuPt30","DiMuon Rapidity after #mu#mu Pt>30 Selection",100,-4,4);
  TH1F* yDiMuPt50 = new TH1F("yDiMuPt50","DiMuon Rapidity after #mu#mu Pt>50 Selection",100,-4,4);
  TH1F* yDiMuPt75 = new TH1F("yDiMuPt75","DiMuon Rapidity after #mu#mu Pt>75 Selection",100,-4,4);

  TH2F* yVptDiMu = new TH2F("yVptDiMu","DiMuon Rapidity v. p_{T}",250,0,500,100,0,4);

  TH1F* ptMu1 = new TH1F("ptMu1","Leading Muon Pt",100,0,200);
  TH1F* ptMu2 = new TH1F("ptMu2","Sub-Leading Muon Pt",100,0,200);
  TH1F* ptJet1 = new TH1F("ptJet1","Leading Jet Pt",200,0,1000);
  TH1F* ptJet2 = new TH1F("ptJet2","Sub-Leading Jet Pt",200,0,1000);

  TH1F* etaMu1 = new TH1F("etaMu1","Leading Muon #eta",50,-2.5,2.5);
  TH1F* etaMu2 = new TH1F("etaMu2","Sub-Leading Muon #eta",50,-2.5,2.5);
  TH1F* etaJet1 = new TH1F("etaJet1","Leading Jet #eta",50,-5.0,5.0);
  TH1F* etaJet2 = new TH1F("etaJet2","Sub-Leading Jet #eta",50,-5.0,5.0);

  TH1F* deltaEtaJets = new TH1F("deltaEtaJets","#Delta#eta Jets",50,0.0,10.0);

  TH1F* countsHist = new TH1F("countsHist","Event Counts",10,0.0,10.0);
  countsHist->GetXaxis()->SetBinLabel(1,"total");
  countsHist->GetXaxis()->SetBinLabel(2,">=2#mu");
  countsHist->GetXaxis()->SetBinLabel(3,"VBFL");
  countsHist->GetXaxis()->SetBinLabel(4,"VBF");
  countsHist->GetXaxis()->SetBinLabel(5,"VBFT");
  countsHist->GetXaxis()->SetBinLabel(6,"Pt30");
  countsHist->GetXaxis()->SetBinLabel(7,"Pt50");
  countsHist->GetXaxis()->SetBinLabel(8,"Pt75");

  TH1F* cosThetaStarHist = new TH1F("cosThetaStar","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarVBFMHist = new TH1F("cosThetaStarVBFM","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarVBFLHist = new TH1F("cosThetaStarVBFL","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarVBFTHist = new TH1F("cosThetaStarVBFT","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarPt30Hist = new TH1F("cosThetaStarPt30","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarPt50Hist = new TH1F("cosThetaStarPt50","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarPt75Hist = new TH1F("cosThetaStarPt75","cos(#theta^{*})",50,-1.,1.);

  TH1F* puJetIDSimpleDiscJet1Hist = new TH1F("puJetIDSimpleDiscJet1","PU Jet ID--Simple Discriminator Leading Jet",50,-1.,1.);
  TH1F* puJetIDSimpleDiscJet2Hist = new TH1F("puJetIDSimpleDiscJet2","PU Jet ID--Simple Discriminator Sub-Leading Jet",50,-1.,1.);
  TH1F* puJetIDSimpleDiscJet3Hist = new TH1F("puJetIDSimpleDiscJet3","PU Jet ID--Simple Discriminator 3rd Leading Jet",50,-1.,1.);

  TH1F* puJetIDSimpleJet1Hist = new TH1F("puJetIDSimpleJet1","PU Jet ID--Simple Loose Leading Jet",2,0,2);
  puJetIDSimpleJet1Hist->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet1Hist->GetXaxis()->SetBinLabel(2,"Pass");
  TH1F* puJetIDSimpleJet2Hist = new TH1F("puJetIDSimpleJet2","PU Jet ID--Simple Loose Sub-Leading Jet",2,-0,2);
  puJetIDSimpleJet2Hist->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet2Hist->GetXaxis()->SetBinLabel(2,"Pass");
  TH1F* puJetIDSimpleJet3Hist = new TH1F("puJetIDSimpleJet3","PU Jet ID--Simple Loose 3rd Leading Jet",2,-0,2);
  puJetIDSimpleJet3Hist->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet3Hist->GetXaxis()->SetBinLabel(2,"Pass");

  TH1F* BDTHistMuonOnly = new TH1F("BDTHistMuonOnly","BDT Discriminator",100,-1,1);
  TH1F* likelihoodHistMuonOnly = new TH1F("likelihoodHistMuonOnly","Likelihood Discriminator",100,-1,1);
  TH1F* LDHistMuonOnly = new TH1F("LDHistMuonOnly","LD Discriminator",100,-1,1);

  TH1F* BDTHistVBF = new TH1F("BDTHistVBF","BDT Discriminator",100,-1,1);
  TH1F* likelihoodHistVBF = new TH1F("likelihoodHistVBF","Likelihood Discriminator",100,-1,1);
  TH1F* LDHistVBF = new TH1F("LDHistVBF","LD Discriminator",100,-1,1);

  //for MVA Applyer

  float cosThetaStar=0.0;
  float mDiJetMVA=-10.0;
  float ptDiJetMVA=-10.0;
  float yDiJetMVA=-10.0;
  float ptMu1MVA=-10.0;
  float ptMu2MVA=-10.0;
  float etaMu1MVA=-10.0;
  float etaMu2MVA=-10.0;
  float ptJet1MVA=-10.0;
  float ptJet2MVA=-10.0;
  float etaJet1MVA=-10.0;
  float etaJet2MVA=-10.0;
  float deltaEtaJetsMVA=-10.0;
  float productEtaJetsMVA=-10.0;
  int nJetsInRapidityGapMVA=-10;

  //Muon Only MVA
   TMVA::Reader * readerMuonOnly = new TMVA::Reader("!Color:!Silent");
   readerMuonOnly->AddVariable( "mDiMu", &recoCandMass);
   readerMuonOnly->AddVariable( "ptDiMu", &recoCandPt);
   readerMuonOnly->AddVariable( "yDiMu", &recoCandY);
   readerMuonOnly->AddSpectator( "mDiJet", &mDiJetMVA);
   readerMuonOnly->AddSpectator( "ptDiJet", &ptDiJetMVA);
   readerMuonOnly->AddSpectator( "yDiJet", &yDiJetMVA);

   readerMuonOnly->AddVariable( "ptMu1", &ptMu1MVA);
   readerMuonOnly->AddVariable( "ptMu2", &ptMu2MVA);
   readerMuonOnly->AddVariable( "etaMu1", &etaMu1MVA);
   readerMuonOnly->AddVariable( "etaMu2", &etaMu2MVA);

   readerMuonOnly->AddSpectator( "ptJet1", &ptJet1MVA);
   readerMuonOnly->AddSpectator( "ptJet2", &ptJet2MVA);
   readerMuonOnly->AddSpectator( "etaJet1", &etaJet1MVA);
   readerMuonOnly->AddSpectator( "etaJet2", &etaJet2MVA);

   readerMuonOnly->AddVariable( "cosThetaStar", &cosThetaStar);
   readerMuonOnly->AddSpectator( "deltaEtaJets", &deltaEtaJetsMVA);
   readerMuonOnly->AddSpectator( "productEtaJets", &productEtaJetsMVA);
   readerMuonOnly->AddSpectator( "nJetsInRapidityGap", &nJetsInRapidityGapMVA);

   readerMuonOnly->BookMVA("BDT","weightsMuonOnly/TMVAClassification_BDT.weights.xml");
   readerMuonOnly->BookMVA("Likelihood","weightsMuonOnly/TMVAClassification_Likelihood.weights.xml");
   readerMuonOnly->BookMVA("LD","weightsMuonOnly/TMVAClassification_LD.weights.xml");

  //VBF MVA
   TMVA::Reader * readerVBF = new TMVA::Reader("!Color:!Silent");
   readerVBF->AddVariable( "mDiMu", &recoCandMass);
   readerVBF->AddVariable( "ptDiMu", &recoCandPt);
   readerVBF->AddVariable( "yDiMu", &recoCandY);
   readerVBF->AddVariable( "mDiJet", &mDiJetMVA);
   readerVBF->AddSpectator( "ptDiJet", &ptDiJetMVA);
   readerVBF->AddSpectator( "yDiJet", &yDiJetMVA);

   readerVBF->AddVariable( "ptMu1", &ptMu1MVA);
   readerVBF->AddVariable( "ptMu2", &ptMu2MVA);
   readerVBF->AddSpectator( "etaMu1", &etaMu1MVA);
   readerVBF->AddSpectator( "etaMu2", &etaMu2MVA);

   readerVBF->AddSpectator( "ptJet1", &ptJet1MVA);
   readerVBF->AddSpectator( "ptJet2", &ptJet2MVA);
   readerVBF->AddSpectator( "etaJet1", &etaJet1MVA);
   readerVBF->AddSpectator( "etaJet2", &etaJet2MVA);

   readerVBF->AddVariable( "cosThetaStar", &cosThetaStar);
   readerVBF->AddVariable( "deltaEtaJets", &deltaEtaJetsMVA);
   readerVBF->AddSpectator( "productEtaJets", &productEtaJetsMVA);
   readerVBF->AddSpectator( "nJetsInRapidityGap", &nJetsInRapidityGapMVA);

   readerVBF->BookMVA("BDT","weightsVBF/TMVAClassification_BDT.weights.xml");
   readerVBF->BookMVA("Likelihood","weightsVBF/TMVAClassification_Likelihood.weights.xml");
   readerVBF->BookMVA("LD","weightsVBF/TMVAClassification_LD.weights.xml");
  
  unsigned nEvents = tree->GetEntries();
  cout << "nEvents: " << nEvents << endl;
  for(unsigned i=0; i<nEvents;i++)
  {
    if(i >= maxEvents)
      continue;
    tree->GetEvent(i);
    if (i % 1000 == 0) cout << "Event: " << i << endl;

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

    cosThetaStar=-10.0;
    mDiJetMVA=-10.0;
    ptDiJetMVA=-10.0;
    yDiJetMVA=-10.0;
    ptMu1MVA=muon1.pt;
    ptMu2MVA=muon2.pt;
    etaMu1MVA=muon1.eta;
    etaMu2MVA=muon2.eta;
    ptJet1MVA=-10.0;
    ptJet2MVA=-10.0;
    etaJet1MVA=-10.0;
    etaJet2MVA=-10.0;
    deltaEtaJetsMVA=-10.0;
    productEtaJetsMVA=-10.0;
    nJetsInRapidityGapMVA=-10;

    //////////////////////////////////////////
    //Computing CosTheta*

    TLorentzVector pMuon1;
    TLorentzVector pMuon2;
    pMuon1.SetPtEtaPhiM(muon1.pt,muon1.eta,muon1.phi,0.105);
    pMuon2.SetPtEtaPhiM(muon2.pt,muon2.eta,muon2.phi,0.105);
    TLorentzVector diMuon = pMuon1+pMuon2;

    TLorentzVector starMuon1 = pMuon1;
    TLorentzVector starMuon2 = pMuon2;
    TVector3 boost = diMuon.BoostVector();
    starMuon1.Boost(-boost);
    starMuon2.Boost(-boost);

    //if (muon1.charge>0)
    if ((int) (1000 * muon1.pt) % 2 == 0)
    {
        TVector3 directionOfBoost = starMuon1.BoostVector();
        cosThetaStar = directionOfBoost.Dot(diMuon.BoostVector()) / (directionOfBoost.Mag()*diMuon.BoostVector().Mag());
    }
    else
    {
        TVector3 directionOfBoost = starMuon2.BoostVector();
        cosThetaStar = directionOfBoost.Dot(diMuon.BoostVector()) / (directionOfBoost.Mag()*diMuon.BoostVector().Mag());
    }

    //////////////////////////////////////////
    // Filling Hists

    mDiMu->Fill(recoCandMass);
    yDiMu->Fill(recoCandY);
    ptDiMu->Fill(recoCandPt);
    yVptDiMu->Fill(recoCandPt,fabs(recoCandY));
    ptMu1->Fill(muon1.pt);
    ptMu2->Fill(muon2.pt);
    etaMu1->Fill(muon1.eta);
    etaMu2->Fill(muon2.eta);
    cosThetaStarHist->Fill(cosThetaStar);
  
    if (recoCandPt>30.0)
    {
      countsHist->Fill(5.0);
      yDiMuPt30->Fill(recoCandY);
      cosThetaStarPt30Hist->Fill(cosThetaStar);
      if (recoCandPt>50.0)
      {
        countsHist->Fill(6.0);
        yDiMuPt50->Fill(recoCandY);
        cosThetaStarPt50Hist->Fill(cosThetaStar);
        if (recoCandPt>75.0)
        {
          countsHist->Fill(7.0);
          yDiMuPt75->Fill(recoCandY);
          cosThetaStarPt75Hist->Fill(cosThetaStar);
        }
      }
    }
    // Eta Categories
    if(fabs(muon1.eta) < 1.0 && fabs(muon2.eta) < 1.0)
    {
      mDiMuEta11->Fill(recoCandMass);
    }
    if(fabs(muon1.eta) < 1.0 || fabs(muon2.eta) < 1.0)
    {
      mDiMuEta12->Fill(recoCandMass);
    }
    else
    {
      mDiMuEta22->Fill(recoCandMass);
    }

    // Jet Part
    bool goodJets = false;
    if(jets.nPFjets>=2 && jets.pfJetPt[0]>30.0 && jets.pfJetPt[1]>30.0)
        goodJets = true;

    bool VBFVeryLoose = false;
    bool vbfLoose = true;
    if(goodJets)
    {
      TLorentzVector pJet1;
      TLorentzVector pJet2;
      pJet1.SetXYZM(jets.pfJetPx[0],jets.pfJetPy[0],jets.pfJetPz[0],jets.pfJetM[0]);
      pJet2.SetXYZM(jets.pfJetPx[1],jets.pfJetPy[1],jets.pfJetPz[1],jets.pfJetM[1]);
      TLorentzVector diJet = pJet1+pJet2;

      mDiJet->Fill(diJet.M());
      double dEtaJets = fabs(jets.pfJetEta[0]-jets.pfJetEta[1]);
      double etaJetProduct = jets.pfJetEta[0]*jets.pfJetEta[1];
      deltaEtaJets->Fill(dEtaJets);
      ptJet1->Fill(jets.pfJetPt[0]);
      ptJet2->Fill(jets.pfJetPt[1]);
      etaJet1->Fill(jets.pfJetEta[0]);
      etaJet2->Fill(jets.pfJetEta[1]);

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
      nJetsInRapidityGapMVA = 0;
      for(unsigned iJet=2; (iJet < jets.nPFjets && iJet < 10);iJet++)
      {
        if(jets.pfJetPt[iJet] > 30.0)
        {
          if(jets.pfJetEta[iJet] < etaMax && jets.pfJetEta[iJet] > etaMin)
          {
            jetInRapidityGap = true;
            nJetsInRapidityGapMVA++;
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

      mDiJetMVA = diJet.M();
      yDiJetMVA = diJet.Rapidity();
      ptDiJetMVA = diJet.Pt();
      ptJet1MVA = pJet1.Pt();
      ptJet2MVA = pJet2.Pt();
      etaJet1MVA = pJet1.Eta();
      etaJet2MVA = pJet2.Eta();
      productEtaJetsMVA = etaJetProduct;

      //VBF MVA
      if (nJetsInRapidityGapMVA==0 && productEtaJetsMVA<0.0)
      {
        VBFVeryLoose=true;
        float likelihooodDisc = readerVBF->EvaluateMVA("Likelihood");
        float BDTDisc = readerVBF->EvaluateMVA("BDT");
        float LDDisc = readerVBF->EvaluateMVA("LD");
        //std::cout << "BDT: "<< BDTDisc << " likelihood: " << likelihooodDisc << " LD: " << LDDisc << std::endl;
        likelihoodHistVBF->Fill(likelihooodDisc);
        LDHistVBF->Fill(LDDisc);
        BDTHistVBF->Fill(BDTDisc);
        mDiMuVBFVL->Fill(recoCandMass);
      }

      //VBFLoose Selection
      vbfLoose = true;
      if(dEtaJets <= 3.0)
          vbfLoose =false;
      if(etaJetProduct >= 0)
          vbfLoose =false;
      if(diJet.M()<=300.0)
          vbfLoose =false;
      if (jetInRapidityGap)
          vbfLoose =false;

      //VBF Selection
      bool vbfSelected = vbfLoose;
      if(dEtaJets <= 4.0)
        vbfSelected=false;
      if(diJet.M()<=400.0)
        vbfSelected=false;
      // Below already checked
      //if(etaJetProduct >= 0)
      //    continue;
      //if (jetInRapidityGap)
      //    continue;
  
      //VBF Tight Selection
      bool vbfTight = vbfSelected;
      if(dEtaJets <= 5.0)
        vbfTight=false;

  
  
//HIG-12-007 PAS H->tautau
//The VBF category requires at least two jets with pT > 30 GeV/c, |η1 − η2 | > 4.0 and
//η1 · η2 < 0 (with η1 the pseudorapidty of the leading jet and η2 the pseudorapidity
//of the subleading jet), and a di-jet invariant mass m12 > 400 GeV/c2 , with no other
//jet with pT > 30 GeV/c in the rapidity region between the two jets.
  
  
      if(vbfTight)
      {
        countsHist->Fill(4.0);
        mDiMuVBFT->Fill(recoCandMass);
        cosThetaStarVBFTHist->Fill(cosThetaStar);
      }
      else if(vbfSelected)
      {
        mDiMuVBFM->Fill(recoCandMass);
        ptDiMuVBFM->Fill(recoCandPt);
        yDiMuVBFM->Fill(recoCandY);
        cosThetaStarVBFMHist->Fill(cosThetaStar);
        countsHist->Fill(3.0);
      }
      else if(vbfLoose)
      {
        mDiMuVBFL->Fill(recoCandMass);
        ptDiMuVBFL->Fill(recoCandPt);
        yDiMuVBFL->Fill(recoCandY);
        cosThetaStarVBFLHist->Fill(cosThetaStar);
        countsHist->Fill(2.0);
      }
  
    } //end jet part

    if (!vbfLoose)
    {
      if(recoCandPt<30.0)
          mDiMuPtL30->Fill(recoCandMass);
      else if(recoCandPt>=30.0 && recoCandPt<50.0)
          mDiMuPt30to50->Fill(recoCandMass);
      else if(recoCandPt>=50.0 && recoCandPt<75.0)
          mDiMuPt50to75->Fill(recoCandMass);
      else if(recoCandPt>=75.0 && recoCandPt<125.0)
          mDiMuPt75to125->Fill(recoCandMass);
      else if(recoCandPt>=125.0)
          mDiMuPt125->Fill(recoCandMass);
    }

    if (!VBFVeryLoose)
    {
      // Muon Only MVA
      if (recoCandMass>110.0 && recoCandMass < 140.0)
      {
        float likelihooodDisc = readerMuonOnly->EvaluateMVA("Likelihood");
        float BDTDisc = readerMuonOnly->EvaluateMVA("BDT");
        float LDDisc = readerMuonOnly->EvaluateMVA("LD");
        //std::cout << "BDT: "<< BDTDisc << " likelihood: " << likelihooodDisc << " LD: " << LDDisc << std::endl;
        likelihoodHistMuonOnly->Fill(likelihooodDisc);
        LDHistMuonOnly->Fill(LDDisc);
        BDTHistMuonOnly->Fill(BDTDisc);
      }

    }


  }// end event loop

  outFile->cd();

  mDiMu->Write();
  mDiJet->Write();

  ptMu1->Write();
  ptMu2->Write();
  ptJet1->Write();
  ptJet2->Write();

  etaMu1->Write();
  etaMu2->Write();
  etaJet1->Write();
  etaJet2->Write();

  deltaEtaJets->Write();
  mDiMuVBFM->Write();
  mDiMuVBFL->Write();
  mDiMuVBFT->Write();
  mDiMuVBFVL->Write();

  mDiMuEta11->Write();
  mDiMuEta12->Write();
  mDiMuEta22->Write();

  ptDiMu->Write();
  ptDiMuVBFM->Write();
  ptDiMuVBFL->Write();

  mDiMuPtL30->Write();
  mDiMuPt30to50->Write();
  mDiMuPt50to75->Write();
  mDiMuPt75to125->Write();
  mDiMuPt125->Write();

  yDiMu->Write();
  yDiMuVBFM->Write();
  yDiMuVBFL->Write();
  yDiMuPt30->Write();
  yDiMuPt50->Write();
  yDiMuPt75->Write();

  yVptDiMu->Write();

  countsHist->Write();

  cosThetaStarHist->Write();
  cosThetaStarVBFMHist->Write();
  cosThetaStarVBFLHist->Write();
  cosThetaStarVBFTHist->Write();
  cosThetaStarPt30Hist->Write();
  cosThetaStarPt50Hist->Write();
  cosThetaStarPt75Hist->Write();

  puJetIDSimpleDiscJet1Hist->Write();
  puJetIDSimpleDiscJet2Hist->Write();
  puJetIDSimpleDiscJet3Hist->Write();

  puJetIDSimpleJet1Hist->Write();
  puJetIDSimpleJet2Hist->Write();
  puJetIDSimpleJet3Hist->Write();

  BDTHistMuonOnly->Write();
  likelihoodHistMuonOnly->Write();
  LDHistMuonOnly->Write();

  BDTHistVBF->Write();
  likelihoodHistVBF->Write();
  LDHistVBF->Write();

  delete readerMuonOnly;
  delete readerVBF;
  cout << "analyzer done." << endl << endl;
  return 0;
}
