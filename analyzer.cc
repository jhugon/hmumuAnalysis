#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>

#include "PhysicsObjectFunctions.h"

double getDiM(TClonesArray*, double maxEta);

int main(int argc, char *argv[])
{
  if (argc != 2) 
  {
	std::cout << "testReader requires 1 argument:" << std::endl
	<< "testReader <infilename.root>" << std::endl;
	return 1;
  }

  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas("c1");
  TFile* f = new TFile(argv[1]);
  TTree * tree = (TTree*) f->Get("tree");
  TH2F* htemp;
  TH1D *htemp_2;

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

  TH1F* muRelIso = new TH1F("muRelIso","Relative Isolation",20,0,2);
  TH1F* elRelIso = new TH1F("elRelIso","Relative Isolation",20,0,2);
  TH1F* phRelIso = new TH1F("phRelIso","Relative Isolation",20,0,2);

  TH1F* mDiMu = new TH1F("mDiMu","DiMuon Mass",50,0,200);
  TH1F* mDiEl = new TH1F("mDiEl","DiElectron Mass",50,0,200);
  TH1F* mDiJet = new TH1F("mDiJet","DiJet Mass",50,0,200);
  TH1F* mDiPh = new TH1F("mDiPh","DiPhoton Mass",50,0,200);
  
  unsigned nLightJets = 0;
  unsigned nBJets = 0;
  unsigned nBTagLightJets = 0;
  unsigned nBTagBJets = 0;
  unsigned nEvents = tree->GetEntries();
  for(unsigned i=0; i<nEvents;i++)
  {
    tree->GetEvent(i);

    unsigned nMuons = muons->GetEntries();
    unsigned nElectrons = electrons->GetEntries();
    unsigned nPhotons = photons->GetEntries();
    unsigned nJets = jets->GetEntries();

    for(unsigned i=0; i<nMuons; i++)
    {
	TParticle* tmpMuon = (TParticle*) muons->At(i);
	muRelIso->Fill(getIso(tmpMuon)/tmpMuon->Pt());
    }
    for(unsigned i=0; i<nElectrons; i++)
    {
	TParticle* tmpEl = (TParticle*) electrons->At(i);
	elRelIso->Fill(getIso(tmpEl)/tmpEl->Pt());
    }
    for(unsigned i=0; i<nPhotons; i++)
    {
	TParticle* tmpPh = (TParticle*) photons->At(i);
	phRelIso->Fill(getIso(tmpPh)/tmpPh->Pt());
    }
    for(unsigned i=0; i<nJets; i++)
    {
	bool isB = getGenB((TParticle*) jets->At(i));
	bool isBTag = getBTag((TParticle*) jets->At(i));
	if(isB)
	{
	  nBJets++;
	  if(isBTag)
	    nBTagBJets++;
	}
	else
	{
	  nLightJets++;
	  if(isBTag)
	    nBTagLightJets++;
	}
    }
	

    if(nMuons>1)
        mDiMu->Fill(getDiM(muons,2.4));
    if(nElectrons>1)
        mDiEl->Fill(getDiM(electrons,2.5));
    if(nPhotons>1)
        mDiPh->Fill(getDiM(photons,2.5));
    if(nJets>1)
        mDiJet->Fill(getDiM(jets,5.0));
/*
    getBTag();
    getGenB();
    getCaloIso();
    getTrkIso();
*/
  }

  muRelIso->SetLineColor(kRed+1);
  elRelIso->SetLineColor(kGreen+1);
  phRelIso->SetLineColor(kYellow+1);
  
  muRelIso->GetXaxis()->SetTitle("Relative Isolation");
  muRelIso->GetYaxis()->SetTitle("Counts");
  muRelIso->Draw();
  elRelIso->Draw("same");
  phRelIso->Draw("same");
  TLegend* isoLegend = new TLegend(0.7,0.70,0.9,0.9);
  isoLegend->SetFillColor(kWhite);
  isoLegend->AddEntry(muRelIso,"#mu Iso","l");
  isoLegend->AddEntry(elRelIso,"e Iso","l");
  isoLegend->AddEntry(phRelIso,"#gamma Iso","l");
  isoLegend->Draw("same");
  c1->SaveAs("iso.png");

  std::cout << "\n##########################\nb-Tags: \n";
  std::cout << "Mistag Rate: " << float(nBTagLightJets)/float(nLightJets)<< std::endl;
  std::cout << "BTag Eff:    " << float(nBTagBJets)/float(nBJets)<< 
							std::endl << std::endl;

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
