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
  if (argc < 2) 
  {
	std::cout << "testReader requires >=1 argument:" << std::endl
	<< "testReader <infilename.root>" << std::endl;
	return 1;
  }

  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas("c1");
  TChain * tree = new TChain("tree");

  for(unsigned iarg=1; iarg<argc;iarg++)
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

  TH1F* mDiMu = new TH1F("mDiMu","DiMuon Mass",50,0,200);
  TH1F* mDiEl = new TH1F("mDiEl","DiElectron Mass",50,0,200);
  TH1F* mDiJet = new TH1F("mDiJet","DiJet Mass",50,0,200);
  TH1F* mDiPh = new TH1F("mDiPh","DiPhoton Mass",50,0,200);

  TH1F* ptMu1 = new TH1F("ptMu1","Leading Muon Pt",50,0,200);
  TH1F* ptMu2 = new TH1F("ptMu2","Sub-Leading Muon Pt",50,0,200);
  TH1F* ptJet1 = new TH1F("ptJet1","Leading Jet Pt",50,0,200);
  TH1F* ptJet2 = new TH1F("ptJet2","Sub-Leading Jet Pt",50,0,200);

  TH1F* deltaEtaJets = new TH1F("deltaEtaJets","#Delta#eta Jets",50,-10.0,10.0);
  
  unsigned nLightJets = 0;
  unsigned nBJets = 0;
  unsigned nBTagLightJets = 0;
  unsigned nBTagBJets = 0;
  unsigned nEvents = tree->GetEntries();
  cout << "nEvents: " << nEvents << endl;
  for(unsigned i=0; i<nEvents;i++)
  {
    tree->GetEvent(i);
    if (i % 1000 == 0) cout << "Event: " << i << endl;

    unsigned nMuons = muons->GetEntries();
    unsigned nElectrons = electrons->GetEntries();
    unsigned nPhotons = photons->GetEntries();
    unsigned nJets = jets->GetEntries();

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
	

    mDiMu->Fill(getDiM(muons,2.4));
    mDiEl->Fill(getDiM(electrons,2.5));
    mDiPh->Fill(getDiM(photons,2.5));
    mDiJet->Fill(getDiM(jets,5.0));
    deltaEtaJets->Fill(getDeltaEta(jets,5.0));
  }

  c1->Clear();
  mDiMu->Draw();
  c1->SaveAs("mDiMu.png");
  c1->Clear();
  mDiEl->Draw();
  c1->SaveAs("mDiEl.png");
  c1->Clear();
  mDiJet->Draw();
  c1->SaveAs("mDiJet.png");
  c1->Clear();
  mDiPh->Draw();
  c1->SaveAs("mDiPh.png");
  ptMu1->Draw();
  c1->SaveAs("ptMu1.png");
  ptMu2->Draw();
  c1->SaveAs("ptMu2.png");
  ptJet1->Draw();
  c1->SaveAs("ptJet1.png");
  ptJet2->Draw();
  c1->SaveAs("ptJet2.png");

  deltaEtaJets->Draw();
  c1->SaveAs("deltaEtaJets.png");

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
