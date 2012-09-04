#include "mva.h"

MVA::MVA(const std::string outFileName)
{
  if(outFileName.size()>0)
  {
    outFile_ = new TFile(outFileName.c_str(),"RECREATE");
    outFile_->cd();
    outTree_ = new TTree("tree","tree");

    outTree_->Branch("mDiMu",&mDiMu,"mDiMu/F");
    outTree_->Branch("ptDiMu",&ptDiMu,"ptDiMu/F");
    outTree_->Branch("yDiMu",&yDiMu,"yDiMu/F");
    outTree_->Branch("mDiJet",&mDiJet,"mDiJet/F");
    outTree_->Branch("ptDiJet",&ptDiJet,"ptDiJet/F");
    outTree_->Branch("yDiJet",&yDiJet,"yDiJet/F");
    outTree_->Branch("ptMu1",&ptMu1,"ptMu1/F");
    outTree_->Branch("ptMu2",&ptMu2,"ptMu2/F");
    outTree_->Branch("etaMu1",&etaMu1,"etaMu1/F");
    outTree_->Branch("etaMu2",&etaMu2,"etaMu2/F");
    outTree_->Branch("ptJet1",&ptJet1,"ptJet1/F");
    outTree_->Branch("ptJet2",&ptJet2,"ptJet2/F");
    outTree_->Branch("etaJet1",&etaJet1,"etaJet1/F");
    outTree_->Branch("etaJet2",&etaJet2,"etaJet2/F");
    outTree_->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
    outTree_->Branch("deltaEtaJets",&deltaEtaJets,"deltaEtaJets/F");
    outTree_->Branch("productEtaJets",&productEtaJets,"productEtaJets/F");
    outTree_->Branch("nJetsInRapidityGap",&nJetsInRapidityGap,"nJetsInRapidityGap/F");
  
    outTree_->Branch("deltaPhiJets",&deltaPhiJets,"deltaPhiJets/F");
    outTree_->Branch("deltaRJets",&deltaRJets,"deltaRJets/F");
    outTree_->Branch("deltaEtaMuons",&deltaEtaMuons,"deltaEtaMuons/F");
    outTree_->Branch("deltaPhiMuons",&deltaPhiMuons,"deltaPhiMuons/F");
    outTree_->Branch("deltaRMuons",&deltaRMuons,"deltaRMuons/F");

    outTree_->Branch("relIsoMu1",&relIsoMu1,"relIsoMu1/F");
    outTree_->Branch("relIsoMu2",&relIsoMu2,"relIsoMu2/F");
    outTree_->Branch("ht",&ht,"ht/F");
    outTree_->Branch("nJets",&nJets,"nJets/F");
    outTree_->Branch("htInRapidityGap",&htInRapidityGap,"htInRapidityGap/F");
  }

  if(outFileName.size()<1)
  {
    reader_ = new TMVA::Reader("!Color:!Silent");
  
    reader_->AddVariable("mDiMu",&mDiMu);
    reader_->AddVariable("ptDiMu",&ptDiMu);
    reader_->AddVariable("yDiMu",&yDiMu);
    reader_->AddVariable("mDiJet",&mDiJet);
    reader_->AddVariable("ptDiJet",&ptDiJet);
    reader_->AddVariable("yDiJet",&yDiJet);
    reader_->AddVariable("ptMu1",&ptMu1);
    reader_->AddVariable("ptMu2",&ptMu2);
    reader_->AddVariable("etaMu1",&etaMu1);
    reader_->AddVariable("etaMu2",&etaMu2);
    reader_->AddVariable("ptJet1",&ptJet1);
    reader_->AddVariable("ptJet2",&ptJet2);
    reader_->AddVariable("etaJet1",&etaJet1);
    reader_->AddVariable("etaJet2",&etaJet2);
    reader_->AddVariable("cosThetaStar",&cosThetaStar);
    reader_->AddVariable("deltaEtaJets",&deltaEtaJets);
    reader_->AddVariable("productEtaJets",&productEtaJets);
    reader_->AddVariable("nJetsInRapidityGap",&nJetsInRapidityGap);
    
    reader_->AddVariable("deltaPhiJets",&deltaPhiJets);
    reader_->AddVariable("deltaRJets",&deltaRJets);
    reader_->AddVariable("deltaEtaMuons",&deltaEtaMuons);
    reader_->AddVariable("deltaPhiMuons",&deltaPhiMuons);
    reader_->AddVariable("deltaRMuons",&deltaRMuons);
  
    reader_->BookMVA("BDT","weightsMuonOnly/TMVAClassification_BDT.weights.xml");
    reader_->BookMVA("Likelihood","weightsMuonOnly/TMVAClassification_Likelihood.weights.xml");
  }
}

MVA::~MVA()
{
  if (outTree_ != NULL)
  {
    outFile_->cd();
    outTree_->Write();
    delete outTree_;
  }
  if (outFile_ != NULL)
  {
    outFile_->Close();
    delete outFile_;
  }
  if (reader_ != NULL)
  {
    delete reader_;
  }
}

void
MVA::getMVA(float& bdtValue, float& lhValue)
{
  if(reader_!=NULL)
  {
    bdtValue = reader_->EvaluateMVA("BDT");
    lhValue = reader_->EvaluateMVA("Likelihood");
  }
  else
  {
    bdtValue = -1000.0;
    lhValue = -1000.0;
  }
}

void
MVA::resetValues()
{
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
  nJetsInRapidityGap=-10.0;

  deltaEtaMuons=-10.0;
  deltaPhiMuons=-10.0;
  deltaRMuons=-10.0;
  deltaPhiJets=-10.0;
  deltaRJets=-10.0;

  relIsoMu1=-10.0;
  relIsoMu2=-10.0;
  ht=0.0;
  nJets=0.0;
  htInRapidityGap=0.0;
}
