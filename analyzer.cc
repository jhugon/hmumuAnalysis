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

#include "PhysicsObjectFunctions.h"

using namespace std;

double getDiM(TClonesArray*, double maxEta);
double getDeltaEta(TClonesArray* particles, double maxEta);

int main(int argc, char *argv[])
{
  if (argc < 3) 
  {
	std::cout << "testReader requires >=2 argument:" << std::endl
	<< "analyzer <outfilename.root> <infilenames.root...>" << std::endl;
	return 1;
  }

  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas("c1");
  TChain * tree = new TChain("tree");

  TFile * outFile = new TFile(argv[1],"RECREATE");

  long maxEvents = 100000000;
  //maxEvents = 20000;

  for(unsigned iarg=2; iarg<argc;iarg++)
  {
    tree->AddFile(argv[iarg]);
  }

  TClonesArray* muons = new TClonesArray("TParticle");
  TClonesArray* electrons = new TClonesArray("TParticle");
  TClonesArray* photons = new TClonesArray("TParticle");
  TClonesArray* jets = new TClonesArray("TParticle");
  double met=0.0;
  double metPhi=0.0;
  double sumPt=0.0;
  double trueMet=0.0;
  double trueMetPhi=0.0;
  double trueSumPt=0.0;

  tree->SetBranchAddress("muons",&muons);
  tree->SetBranchAddress("electrons",&electrons);
  tree->SetBranchAddress("photons",&photons);
  tree->SetBranchAddress("jets",&jets);
  tree->SetBranchAddress("met",&met);
  tree->SetBranchAddress("metPhi",&metPhi);
  tree->SetBranchAddress("sumPt",&sumPt);
  tree->SetBranchAddress("trueMet",&trueMet);
  tree->SetBranchAddress("trueMetPhi",&trueMetPhi);
  tree->SetBranchAddress("trueSumPt",&trueSumPt);

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

    unsigned nMuons = muons->GetEntries();
    unsigned nElectrons = electrons->GetEntries();
    unsigned nPhotons = photons->GetEntries();
    unsigned nJets = jets->GetEntries();

    countsHist->Fill(0.0);

    if(nMuons <2)
	continue;

    TParticle * muon1 = (TParticle*) muons->At(0);
    TParticle * muon2 = (TParticle*) muons->At(1);
    if(fabs(muon1->Eta())>2.1 || fabs(muon2->Eta())>2.1 || muon1->Pt()<45.0 || muon2->Pt()<45.0 )
         continue;

    countsHist->Fill(1.0);
/*
    for(unsigned i=0; i<nMuons; i++)
    {
	TParticle* tmpMuon = (TParticle*) muons->At(i);
	float relIso = getIso(tmpMuon)/tmpMuon->Pt();
    }
    for(unsigned i=0; i<nElectrons; i++)
    {
	TParticle* tmpEl = (TParticle*) electrons->At(i);
	float relIso = getIso(tmpEl)/tmpEl->Pt();
    }
    for(unsigned i=0; i<nPhotons; i++)
    {
	TParticle* tmpPh = (TParticle*) photons->At(i);
	float relIso = getIso(tmpPh)/tmpPh->Pt();
    }
    for(unsigned i=0; i<nJets; i++)
    {
	bool isB = getGenB((TParticle*) jets->At(i));
	bool isBTag = getBTag((TParticle*) jets->At(i));
    }
*/

    if(nMuons>=2)
    {
      TLorentzVector pMuon1;
      TLorentzVector pMuon2;
      muon1->Momentum(pMuon1);
      muon2->Momentum(pMuon2);
      TLorentzVector diMuon = pMuon1+pMuon2;
      mDiMu->Fill(diMuon.M());
      yDiMu->Fill(diMuon.Rapidity());
      ptDiMu->Fill(diMuon.Pt());
      ptMu1->Fill(muon1->Pt());
      ptMu2->Fill(muon2->Pt());
      etaMu1->Fill(muon1->Eta());
      etaMu2->Fill(muon2->Eta());

      if (diMuon.Pt()>30.0)
      {
        countsHist->Fill(5.0);
        mDiMuZPt30Selected->Fill(diMuon.M());
        yDiMuZPt30Selected->Fill(diMuon.Rapidity());
        if (diMuon.Pt()>50.0)
        {
          countsHist->Fill(6.0);
          mDiMuZPt50Selected->Fill(diMuon.M());
          yDiMuZPt50Selected->Fill(diMuon.Rapidity());
          if (diMuon.Pt()>75.0)
          {
            countsHist->Fill(7.0);
            mDiMuZPt75Selected->Fill(diMuon.M());
            yDiMuZPt75Selected->Fill(diMuon.Rapidity());
          }
        }
      }
    }
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
    TLorentzVector pMuon1;
    TLorentzVector pMuon2;
    muon1->Momentum(pMuon1);
    muon2->Momentum(pMuon2);
    TLorentzVector diMuon = pMuon1+pMuon2;
    mDiMuVBFLooseSelected->Fill(diMuon.M());
    ptDiMuVBFLooseSelected->Fill(diMuon.Pt());
    yDiMuVBFLooseSelected->Fill(diMuon.Rapidity());
    countsHist->Fill(2.0);

/*
HIG-12-007 PAS H->tautau
The VBF category requires at least two jets with pT > 30 GeV/c, |η1 − η2 | > 4.0 and
η1 · η2 < 0 (with η1 the pseudorapidty of the leading jet and η2 the pseudorapidity
of the subleading jet), and a di-jet invariant mass m12 > 400 GeV/c2 , with no other
jet with pT > 30 GeV/c in the rapidity region between the two jets.
*/


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
    mDiMuVBFSelected->Fill(diMuon.M());
    ptDiMuVBFSelected->Fill(diMuon.Pt());
    yDiMuVBFSelected->Fill(diMuon.Rapidity());
    countsHist->Fill(3.0);

    //VBF Tight Selection
    if(fabs(jet1->Eta()-jet2->Eta())<= 5.0)
        continue;
    //if(diJet.M()<=500.0)
    //    continue;
    countsHist->Fill(4.0);
    mDiMuVBFTightSelected->Fill(diMuon.M());

  }

  c1->Clear();
  mDiMu->Draw();
  c1->SaveAs("mDiMu.png");
  ptDiMu->Draw();
  c1->SaveAs("ptDiMu.png");
  mDiJet->Draw();
  c1->SaveAs("mDiJet.png");

  ptMu1->Draw();
  c1->SaveAs("ptMu1.png");
  ptMu2->Draw();
  c1->SaveAs("ptMu2.png");
  ptJet1->Draw();
  c1->SaveAs("ptJet1.png");
  ptJet2->Draw();
  c1->SaveAs("ptJet2.png");

  etaMu1->Draw();
  c1->SaveAs("etaMu1.png");
  etaMu2->Draw();
  c1->SaveAs("etaMu2.png");
  etaJet1->Draw();
  c1->SaveAs("etaJet1.png");
  etaJet2->Draw();
  c1->SaveAs("etaJet2.png");

  deltaEtaJets->Draw();
  c1->SaveAs("deltaEtaJets.png");

  mDiMuVBFSelected->Draw();
  c1->SaveAs("mDiMuVBFSelected.png");
  ptDiMuVBFSelected->Draw();
  c1->SaveAs("ptDiMuVBFSelected.png");

  mDiMuVBFLooseSelected->Draw();
  c1->SaveAs("mDiMuVBFLooseSelected.png");
  mDiMuVBFTightSelected->Draw();
  c1->SaveAs("mDiMuVBFTightSelected.png");
  ptDiMuVBFLooseSelected->Draw();
  c1->SaveAs("ptDiMuVBFLooseSelected.png");

  mDiMuZPt30Selected->Draw();
  c1->SaveAs("mDiMuZPt30Selected.png");
  mDiMuZPt50Selected->Draw();
  c1->SaveAs("mDiMuZPt50Selected.png");
  mDiMuZPt75Selected->Draw();
  c1->SaveAs("mDiMuZPt75Selected.png");

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

  return 0;
}

double getDiM(TClonesArray* particles, double maxEta)
{
        if(particles->GetEntries() < 2)
            return -1.0;
        TParticle* part0 = NULL;
        TParticle* part1 = NULL;
        for(unsigned i=0; i< particles->GetEntries();i++)
        {
            TParticle * tmpPart = (TParticle*) particles->At(0);
            if(fabs(tmpPart->Eta())<maxEta)
            {
                if(part0 == NULL)
                {
                    part0 = tmpPart;
                }
                else if(part1 == NULL)
                {
                    part1 = tmpPart;
                    break;
                }
            }
        }
        if(part0 == NULL || part1 == NULL)
            return -1.0;

        TLorentzVector p0;
        TLorentzVector p1;
        part0->Momentum(p0);
        part1->Momentum(p1);
        return (p0+p1).M();
}

double getDeltaEta(TClonesArray* particles, double maxEta)
{
        if(particles->GetEntries() < 2)
            return -999999.0;
        TParticle* part0 = NULL;
        TParticle* part1 = NULL;
        for(unsigned i=0; i< particles->GetEntries();i++)
        {
            TParticle * tmpPart = (TParticle*) particles->At(0);
            if(fabs(tmpPart->Eta())<maxEta)
            {
                if(part0 == NULL)
                {
                    part0 = tmpPart;
                }
                else if(part1 == NULL)
                {
                    part1 = tmpPart;
                    break;
                }
            }
        }
        if(part0 == NULL || part1 == NULL)
            return -999999.0;
        return part0->Eta() - part1->Eta();
}
