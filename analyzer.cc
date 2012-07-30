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

  TFile * outFile = new TFile("outfile.root","RECREATE");

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
  TH1F* mDiMuSelected = new TH1F("mDiMuSelected","DiMuon Mass after VBF Selection",50,0,200);
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
    TParticle * jet1 = (TParticle*) jets->At(0);
    TParticle * jet2 = (TParticle*) jets->At(1);
    if(fabs(jet1->Eta()-jet2->Eta())<= 4.0)
        continue;
    if(jet1->Eta()*jet2->Eta() >= 0)
        continue;

    TLorentzVector pJet1;
    TLorentzVector pJet2;
    jet1->Momentum(pJet1);
    jet2->Momentum(pJet2);
    TLorentzVector diJet = pJet1+pJet2;
    if(diJet.M()<=400.0)
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
    //VBF Selected
    mDiMuSelected->Fill(getDiM(muons,2.4));
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

  mDiMuSelected->Draw();
  c1->SaveAs("mDiMuSelected.png");

  outFile->cd();

  mDiMu->Write();
  mDiEl->Write();
  mDiJet->Write();
  mDiPh->Write();
  ptMu1->Write();
  ptMu2->Write();
  ptJet1->Write();
  ptJet2->Write();
  deltaEtaJets->Write();
  mDiMuSelected->Write();

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
