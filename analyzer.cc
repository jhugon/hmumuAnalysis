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
	std::cout << "testReader requires >=2 argument:" << std::endl
	<< "analyzer <outfilename.root> <infilenames.root...>" << std::endl;
	return 1;
  }


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

  //////////////////////////
  // Histograms

  TH1F* mDiMu = new TH1F("mDiMu","DiMuon Mass",200,0,200);
  TH1F* mDiMuVBFSelected = new TH1F("mDiMuVBFSelected","DiMuon Mass after VBF Selection",200,0,200);
  TH1F* mDiMuVBFLooseSelected = new TH1F("mDiMuVBFLooseSelected","DiMuon Mass after VBFLoose Selection",200,0,200);
  TH1F* mDiMuVBFTightSelected = new TH1F("mDiMuVBFTightSelected","DiMuon Mass after VBFTight Selection",200,0,200);
  TH1F* mDiMuZPt30Selected = new TH1F("mDiMuZPt30Selected","DiMuon Mass after p_T^{#mu#mu}>30 GeV Selection",200,0,200);
  TH1F* mDiMuZPt50Selected = new TH1F("mDiMuZPt50Selected","DiMuon Mass after p_T^{#mu#mu}>50 GeV Selection",200,0,200);
  TH1F* mDiMuZPt75Selected = new TH1F("mDiMuZPt75Selected","DiMuon Mass after p_T^{#mu#mu}>75 GeV Selection",200,0,200);
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

    countsHist->Fill(1.0);

    //////////////////////////////////////////
    //Computing CosTheta*

    TLorentzVector pMuon1;
    TLorentzVector pMuon2;
    pMuon1.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,0.105);
    pMuon2.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,0.105);
    TLorentzVector diMuon = pMuon1+pMuon2;

    TLorentzVector starMuon1 = pMuon1;
    TLorentzVector starMuon2 = pMuon2;
    TVector3 boost = diMuon.BoostVector();
    starMuon1.Boost(-boost);
    starMuon2.Boost(-boost);


    double cosThetaStar=0.0;
    if (reco1.charge>0)
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
    ptMu1->Fill(reco1.pt);
    ptMu2->Fill(reco2.pt);
    etaMu1->Fill(reco1.eta);
    etaMu2->Fill(reco2.eta);
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
/*
    if(nJets>=2)
    {
      TParticle * jet1 = (TParticle*) jets->At(0);
      TParticle * jet2 = (TParticle*) jets->At(1);
  
      TLorentzVector pJet1;
      TLorentzVector pJet2;
      jet1->Momentum(pJet1);
      jet2->Momentum(pJet2);
      TLorentzVector diJet = pJet1+pJet2;

      mDiJet->Fill(diJet.M());
      deltaEtaJets->Fill(fabs(jet1->Eta()-jet2->Eta()));
      ptJet1->Fill(jet1->Pt());
      ptJet2->Fill(jet2->Pt());
      etaJet1->Fill(jet1->Eta());
      etaJet2->Fill(jet2->Eta());
    }

    //VBFLoose Selection
    if(muons->GetEntries()<2)
        continue;
    if(jets->GetEntries()<2)
        continue;
    TParticle * jet1 = (TParticle*) jets->At(0);
    TParticle * jet2 = (TParticle*) jets->At(1);
    if(fabs(jet1->Eta()-jet2->Eta())<= 3.0)
        continue;
    if(jet1->Eta()*jet2->Eta() >= 0)
        continue;

    TLorentzVector pJet1;
    TLorentzVector pJet2;
    jet1->Momentum(pJet1);
    jet2->Momentum(pJet2);
    TLorentzVector diJet = pJet1+pJet2;
    if(diJet.M()<=300.0)
        continue;

    float etaMax = jet1->Eta();
    float etaMin = 9999999.0;
    if(etaMax < jet2->Eta())
    {
        etaMax = jet2->Eta();
        etaMin = jet1->Eta();
    }
    else
    {
        etaMin = jet2->Eta();
    }
    bool jetInRapidityGap=false;
    for(float iJet=2; iJet<jets->GetEntries();iJet++)
    {
      TParticle * jet = (TParticle*) jets->At(0);
      if(jet->Eta() < etaMax && jet->Eta() > etaMin)
      {
        jetInRapidityGap = true;
        break;
      }
    }
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
    if(muons->GetEntries()<2)
        continue;
    if(jets->GetEntries()<2)
        continue;
    //TParticle * jet1 = (TParticle*) jets->At(0);
    //TParticle * jet2 = (TParticle*) jets->At(1);
    if(fabs(jet1->Eta()-jet2->Eta())<= 4.0)
        continue;
    if(jet1->Eta()*jet2->Eta() >= 0)
        continue;

    //TLorentzVector pJet1;
    //TLorentzVector pJet2;
    //jet1->Momentum(pJet1);
    //jet2->Momentum(pJet2);
    //TLorentzVector diJet = pJet1+pJet2;
    if(diJet.M()<=400.0)
        continue;

    //float etaMax = jet1->Eta();
    //float etaMin = 9999999.0;
    //if(etaMax < jet2->Eta())
    //{
    //    etaMax = jet2->Eta();
    //    etaMin = jet1->Eta();
    //}
    //else
    //{
    //    etaMin = jet2->Eta();
    //}
    //bool jetInRapidityGap=false;
    //for(float iJet=2; iJet<jets->GetEntries();iJet++)
    //{
    //  TParticle * jet = (TParticle*) jets->At(0);
    //  if(jet->Eta() < etaMax && jet->Eta() > etaMin)
    //  {
    //    jetInRapidityGap = true;
    //    break;
    //  }
    //}
    //if (jetInRapidityGap)
    //    continue;
    //VBF Selected
    mDiMuVBFSelected->Fill(recoCandMass);
    ptDiMuVBFSelected->Fill(recoCandPt);
    yDiMuVBFSelected->Fill(recoCandY);
    cosThetaStarVBFSelectedHist->Fill(cosThetaStar);
    countsHist->Fill(3.0);

    //VBF Tight Selection
    if(fabs(jet1->Eta()-jet2->Eta())<= 5.0)
        continue;
    //if(diJet.M()<=500.0)
    //    continue;
    countsHist->Fill(4.0);
    mDiMuVBFTightSelected->Fill(recoCandMass);
    cosThetaStarVBFTightSelectedHist->Fill(cosThetaStar);

*/
  }

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

  countsHist->Write();

  cosThetaStarHist->Write();
  cosThetaStarVBFSelectedHist->Write();
  cosThetaStarVBFLooseSelectedHist->Write();
  cosThetaStarVBFTightSelectedHist->Write();
  cosThetaStarZPt30SelectedHist->Write();
  cosThetaStarZPt50SelectedHist->Write();
  cosThetaStarZPt75SelectedHist->Write();

  cout << "analyzer done." << endl << endl;
  return 0;
}
