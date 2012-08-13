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
//  minMmm = 123.0;
//  maxMmm = 127.0;

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

  TH1F* mDiMu = new TH1F("mDiMu","DiMuon Mass",400,0,400);
  TH1F* mDiMuVBFSelected = new TH1F("mDiMuVBFSelected","DiMuon Mass after VBF Selection",400,0,400);
  TH1F* mDiMuVBFLooseSelected = new TH1F("mDiMuVBFLooseSelected","DiMuon Mass after VBFLoose Selection",400,0,400);
  TH1F* mDiMuVBFTightSelected = new TH1F("mDiMuVBFTightSelected","DiMuon Mass after VBFTight Selection",400,0,400);
  TH1F* mDiMuZPt30Selected = new TH1F("mDiMuZPt30Selected","DiMuon Mass after p_T^{#mu#mu}>30 GeV Selection",400,0,400);
  TH1F* mDiMuZPt50Selected = new TH1F("mDiMuZPt50Selected","DiMuon Mass after p_T^{#mu#mu}>50 GeV Selection",400,0,400);
  TH1F* mDiMuZPt75Selected = new TH1F("mDiMuZPt75Selected","DiMuon Mass after p_T^{#mu#mu}>75 GeV Selection",400,0,400);

  TH1F* mDiMuEta11 = new TH1F("mDiMuEta11","DiMuon Mass",400,0,400);
  TH1F* mDiMuEta12 = new TH1F("mDiMuEta12","DiMuon Mass",400,0,400);
  TH1F* mDiMuEta22 = new TH1F("mDiMuEta22","DiMuon Mass",400,0,400);

  TH1F* mDiJet = new TH1F("mDiJet","DiJet Mass",500,0,2000);

  TH1F* ptDiMu = new TH1F("ptDiMu","DiMuon Pt",250,0,500);
  TH1F* ptDiMuVBFSelected = new TH1F("ptDiMuVBFSelected","DiMuon Pt after VBF Selection",250,0,500);
  TH1F* ptDiMuVBFLooseSelected = new TH1F("ptDiMuVBFLooseSelected","DiMuon Pt after VBFLoose Selection",250,0,500);

  TH1F* yDiMu = new TH1F("yDiMu","DiMuon Rapidity",100,-4,4);
  TH1F* yDiMuVBFSelected = new TH1F("yDiMuVBFSelected","DiMuon Rapidity after VBF Selection",100,-4,4);
  TH1F* yDiMuVBFLooseSelected = new TH1F("yDiMuVBFLooseSelected","DiMuon Rapidity after VBFLoose Selection",100,-4,4);
  TH1F* yDiMuZPt30Selected = new TH1F("yDiMuZPt30Selected","DiMuon Rapidity after #mu#mu Pt>30 Selection",100,-4,4);
  TH1F* yDiMuZPt50Selected = new TH1F("yDiMuZPt50Selected","DiMuon Rapidity after #mu#mu Pt>50 Selection",100,-4,4);
  TH1F* yDiMuZPt75Selected = new TH1F("yDiMuZPt75Selected","DiMuon Rapidity after #mu#mu Pt>75 Selection",100,-4,4);

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
  countsHist->GetXaxis()->SetBinLabel(6,"ZPt30");
  countsHist->GetXaxis()->SetBinLabel(7,"ZPt50");
  countsHist->GetXaxis()->SetBinLabel(8,"ZPt75");

  TH1F* cosThetaStarHist = new TH1F("cosThetaStar","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarVBFSelectedHist = new TH1F("cosThetaStarVBFSelected","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarVBFLooseSelectedHist = new TH1F("cosThetaStarVBFLooseSelected","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarVBFTightSelectedHist = new TH1F("cosThetaStarVBFTightSelected","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarZPt30SelectedHist = new TH1F("cosThetaStarZPt30Selected","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarZPt50SelectedHist = new TH1F("cosThetaStarZPt50Selected","cos(#theta^{*})",50,-1.,1.);
  TH1F* cosThetaStarZPt75SelectedHist = new TH1F("cosThetaStarZPt75Selected","cos(#theta^{*})",50,-1.,1.);

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
  
  unsigned nLightJets = 0;
  unsigned nBJets = 0;
  unsigned nBTagLightJets = 0;
  unsigned nBTagBJets = 0;
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


    double cosThetaStar=0.0;
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
      mDiMuZPt30Selected->Fill(recoCandMass);
      yDiMuZPt30Selected->Fill(recoCandY);
      cosThetaStarZPt30SelectedHist->Fill(cosThetaStar);
      if (recoCandPt>50.0)
      {
        countsHist->Fill(6.0);
        mDiMuZPt50Selected->Fill(recoCandMass);
        yDiMuZPt50Selected->Fill(recoCandY);
        cosThetaStarZPt50SelectedHist->Fill(cosThetaStar);
        if (recoCandPt>75.0)
        {
          countsHist->Fill(7.0);
          mDiMuZPt75Selected->Fill(recoCandMass);
          yDiMuZPt75Selected->Fill(recoCandY);
          cosThetaStarZPt75SelectedHist->Fill(cosThetaStar);
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
      for(unsigned iJet=2; (iJet < jets.nPFjets && iJet < 10);iJet++)
      {
        if(jets.pfJetPt[iJet] > 30.0)
        {
          if(jets.pfJetEta[iJet] < etaMax && jets.pfJetEta[iJet] > etaMin)
          {
            jetInRapidityGap = true;
            break;
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

      //VBFLoose Selection
      if(dEtaJets <= 3.0)
          continue;
      if(etaJetProduct >= 0)
          continue;
      if(diJet.M()<=300.0)
          continue;
      if (jetInRapidityGap)
          continue;
  
      //VBFLoose Selected
      mDiMuVBFLooseSelected->Fill(recoCandMass);
      ptDiMuVBFLooseSelected->Fill(recoCandPt);
      yDiMuVBFLooseSelected->Fill(recoCandY);
      cosThetaStarVBFLooseSelectedHist->Fill(cosThetaStar);
      countsHist->Fill(2.0);
  
//HIG-12-007 PAS H->tautau
//The VBF category requires at least two jets with pT > 30 GeV/c, |η1 − η2 | > 4.0 and
//η1 · η2 < 0 (with η1 the pseudorapidty of the leading jet and η2 the pseudorapidity
//of the subleading jet), and a di-jet invariant mass m12 > 400 GeV/c2 , with no other
//jet with pT > 30 GeV/c in the rapidity region between the two jets.
  
  
      //VBF Selection
      if(dEtaJets <= 4.0)
          continue;
      if(diJet.M()<=400.0)
          continue;
      // Below already checked
      //if(etaJetProduct >= 0)
      //    continue;
      //if (jetInRapidityGap)
      //    continue;
  
      //VBF Selected
      mDiMuVBFSelected->Fill(recoCandMass);
      ptDiMuVBFSelected->Fill(recoCandPt);
      yDiMuVBFSelected->Fill(recoCandY);
      cosThetaStarVBFSelectedHist->Fill(cosThetaStar);
      countsHist->Fill(3.0);
  
      //VBF Tight Selection
      if(dEtaJets <= 5.0)
          continue;

      //VBF Tight Selected
      countsHist->Fill(4.0);
      mDiMuVBFTightSelected->Fill(recoCandMass);
      cosThetaStarVBFTightSelectedHist->Fill(cosThetaStar);
    } //end jet part
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
  mDiMuVBFSelected->Write();
  mDiMuVBFLooseSelected->Write();
  mDiMuVBFTightSelected->Write();

  mDiMuEta11->Write();
  mDiMuEta12->Write();
  mDiMuEta22->Write();

  ptDiMu->Write();
  ptDiMuVBFSelected->Write();
  ptDiMuVBFLooseSelected->Write();

  mDiMuZPt30Selected->Write();
  mDiMuZPt50Selected->Write();
  mDiMuZPt75Selected->Write();

  yDiMu->Write();
  yDiMuVBFSelected->Write();
  yDiMuVBFLooseSelected->Write();
  yDiMuZPt30Selected->Write();
  yDiMuZPt50Selected->Write();
  yDiMuZPt75Selected->Write();

  yVptDiMu->Write();

  countsHist->Write();

  cosThetaStarHist->Write();
  cosThetaStarVBFSelectedHist->Write();
  cosThetaStarVBFLooseSelectedHist->Write();
  cosThetaStarVBFTightSelectedHist->Write();
  cosThetaStarZPt30SelectedHist->Write();
  cosThetaStarZPt50SelectedHist->Write();
  cosThetaStarZPt75SelectedHist->Write();

  puJetIDSimpleDiscJet1Hist->Write();
  puJetIDSimpleDiscJet2Hist->Write();
  puJetIDSimpleDiscJet3Hist->Write();

  puJetIDSimpleJet1Hist->Write();
  puJetIDSimpleJet2Hist->Write();
  puJetIDSimpleJet3Hist->Write();

  cout << "analyzer done." << endl << endl;
  return 0;
}
