#ifndef hmumu_mva_h
#define hmumu_mva_h

#include <map>

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
#include "helpers.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace std;

class MVA
{
  public:


  MVA(const std::vector<std::string> configFileNames, const std::string outFileName);
  ~MVA();
  float getMVA(const std::string configFileName, const std::string mvaName);
  void resetValues();

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
  float nJetsInRapidityGap;

  float deltaEtaMuons;
  float deltaPhiMuons;
  float deltaRMuons;
  float deltaPhiJets;
  float deltaRJets;

  float relIsoMu1;
  float relIsoMu2;
  float ht;
  float nJets;
  float htInRapidityGap;

  //////////////////////
  
  std::map<std::string,TMVA::Reader*> readers_;
  TFile* outFile_;
  TTree* outTree_;

  void writeEvent(){if(outTree_ != NULL) outTree_->Fill();};

};

#endif
