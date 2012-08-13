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

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 3) 
  {
	std::cout << "trainingTreeMaker requires >=2 argument:" << std::endl
	<< "trainingTreeMaker <outfilename.root> <infilenames.root...>" << std::endl;
	return 1;
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1");

  TChain * tree = new TChain("tree");

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

  //////////////////////////////////////////////////

  float mDiMu=-10.0;
  float ptDiMu=-10.0;
  float yDiMu=-10.0;
  float mDiJet=-10.0;
  float ptDiJet=-10.0;
  float yDiJet=-10.0;
  float ptMu1=-10.0;
  float ptMu2=-10.0;
  float etaMu1=-10.0;
  float etaMu2=-10.0;
  float ptJet1=-10.0;
  float ptJet2=-10.0;
  float etaJet1=-10.0;
  float etaJet2=-10.0;
  float cosThetaStar=-10.0;
  float deltaEtaJets=-10.0;
  float productEtaJets=-10.0;
  int nJetsInRapidityGap=-10;

  TFile * outFile = new TFile(argv[1],"RECREATE");
  outFile->cd();

  TTree * outTree = new TTree("tree","tree");

  outTree->Branch("mDiMu",&mDiMu,"mDiMu/F");
  outTree->Branch("ptDiMu",&ptDiMu,"ptDiMu/F");
  outTree->Branch("yDiMu",&yDiMu,"yDiMu/F");
  outTree->Branch("mDiJet",&mDiJet,"mDiJet/F");
  outTree->Branch("ptDiJet",&ptDiJet,"ptDiJet/F");
  outTree->Branch("yDiJet",&yDiJet,"yDiJet/F");
  outTree->Branch("ptMu1",&ptMu1,"ptMu1/F");
  outTree->Branch("ptMu2",&ptMu2,"ptMu2/F");
  outTree->Branch("etaMu1",&etaMu1,"etaMu1/F");
  outTree->Branch("etaMu2",&etaMu2,"etaMu2/F");
  outTree->Branch("ptJet1",&ptJet1,"ptJet1/F");
  outTree->Branch("ptJet2",&ptJet2,"ptJet2/F");
  outTree->Branch("etaJet1",&etaJet1,"etaJet1/F");
  outTree->Branch("etaJet2",&etaJet2,"etaJet2/F");
  outTree->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
  outTree->Branch("deltaEtaJets",&deltaEtaJets,"deltaEtaJets/F");
  outTree->Branch("productEtaJets",&productEtaJets,"productEtaJets/F");
  outTree->Branch("nJetsInRapidityGap",&nJetsInRapidityGap,"nJetsInRapidityGap/I");
  
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

    mDiMu=-10.0;
    ptDiMu=-10.0;
    yDiMu=-10.0;
    mDiJet=-10.0;
    ptDiJet=-10.0;
    yDiJet=-10.0;
    ptMu1=-10.0;
    ptMu2=-10.0;
    etaMu1=-10.0;
    etaMu2=-10.0;
    ptJet1=-10.0;
    ptJet2=-10.0;
    etaJet1=-10.0;
    etaJet2=-10.0;
    cosThetaStar=-10.0;
    deltaEtaJets=-10.0;
    productEtaJets=-10.0;
    nJetsInRapidityGap=-10;

    if (!isKinTight_2012(reco1) || !isKinTight_2012(reco2))
        continue;

    //if (recoCandMass < minMmm || recoCandMass > maxMmm)
    //    continue;

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

    mDiMu = recoCandMass;
    yDiMu = recoCandY;
    ptDiMu = recoCandPt;
    ptMu1 = muon1.pt;
    ptMu2 = muon2.pt;
    etaMu1 = muon1.eta;
    etaMu2 = muon2.eta;
  
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

      mDiJet = diJet.M();
      ptDiJet = diJet.Pt();
      yDiJet = diJet.Rapidity();
      deltaEtaJets = fabs(jets.pfJetEta[0]-jets.pfJetEta[1]);
      productEtaJets = jets.pfJetEta[0]*jets.pfJetEta[1];
      ptJet1 = jets.pfJetPt[0];
      ptJet2 = jets.pfJetPt[1];
      etaJet1 = jets.pfJetEta[0];
      etaJet2 = jets.pfJetEta[1];

      // Seeing if there are jets in the rapidity gap
      nJetsInRapidityGap = 0;
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
            nJetsInRapidityGap +=1;
          }
        }
      }

    } //end jet part


    outTree->Fill();

  }

  outFile->cd();
  outTree->Write();

  delete outTree;
  outFile->Close();
  delete tree;

  return 0;
}
