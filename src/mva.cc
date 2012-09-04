#include "mva.h"

MVA::MVA(std::string outFileName)
{
  if(outFileName.size()>0)
  {
    outFile_ = new TFile(outFileName.c_str(),"RECREATE");
    outFile_->cd();
    outTree_ = new TTree("tree","tree");

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
  
    outTree->Branch("deltaPhiJets",&deltaPhiJets,"deltaPhiJets/F");
    outTree->Branch("deltaRJets",&deltaRJets,"deltaRJets/F");
    outTree->Branch("deltaEtaMuons",&deltaEtaMuons,"deltaEtaMuons/F");
    outTree->Branch("deltaPhiMuons",&deltaPhiMuons,"deltaPhiMuons/F");
    outTree->Branch("deltaRMuons",&deltaRMuons,"deltaRMuons/F");
  }

  reader_ = new TMVA::Reader("!Color:!Silent");
}

MVA::~MVA()
{

}

MVA::getMVA()
{

}
