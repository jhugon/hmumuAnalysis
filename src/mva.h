#ifndef hmumu_mva_h
#define hmumu_mva_h

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

class MVA
{
  public:

  MVA();
  ~MVA();
  float getMVA();

  /////////////////////

  float mDiMu;
  float ptDiMu;
  float yDiMu;
  float mDiJet;
  float ptDiJet;
  float yDiJet;
  float ptMu1;
  float ptMu2;
  float etaMu1;
  float etaMu2;
  float ptJet1;
  float ptJet2;
  float etaJet1;
  float etaJet2;
  float cosThetaStar;
  float deltaEtaJets;
  float productEtaJets;
  int nJetsInRapidityGap;

  float deltaEtaMuons;
  float deltaPhiMuons;
  float deltaRMuons;
  float deltaPhiJets;
  float deltaRJets;

  //////////////////////
  
  TMVA::Reader * reader_;
  TTree * outTree_;
  TTree * outFile_;
}

#endif
