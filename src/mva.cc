#include "mva.h"
#include "boost/program_options.hpp"
#include <fstream>

MVA::MVA(const std::vector<std::string> configFileNames, const std::string outFileName)
{
  outTree_ = NULL;
  outFile_ = NULL;

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
    outTree_->Branch("cosThetaStarCS",&cosThetaStarCS,"cosThetaStarCS/F");
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

    outTree_->Branch("nVtx",&nVtx,"nVtx/F");

    outTree_->Branch("weight",&weight,"weight/F");
  }

  if(outFileName.size()<1)
  {
    using namespace boost;

    std::vector<std::string>::const_iterator configFileName = configFileNames.begin();
    std::vector<std::string>::const_iterator configFileNameEnd = configFileNames.end();
    for(;configFileName != configFileNameEnd; configFileName++)
    {
      std::ifstream infile(configFileName->c_str());
      if(!infile)
      {
        std::cout << "Error: Can't open config file...exiting." << std::endl;
        throw;
      }
  
      program_options::options_description optionDesc;
      optionDesc.add_options()
          ("mDiMu",program_options::value<int>(),"")
          ("ptDiMu",program_options::value<int>(),"")
          ("yDiMu",program_options::value<int>(),"")
          ("mDiJet",program_options::value<int>(),"")
          ("ptDiJet",program_options::value<int>(),"")
          ("yDiJet",program_options::value<int>(),"")
          ("ptMu1",program_options::value<int>(),"")
          ("ptMu2",program_options::value<int>(),"")
          ("etaMu1",program_options::value<int>(),"")
          ("etaMu2",program_options::value<int>(),"")
          ("ptJet1",program_options::value<int>(),"")
          ("ptJet2",program_options::value<int>(),"")
          ("etaJet1",program_options::value<int>(),"")
          ("etaJet2",program_options::value<int>(),"")
          ("cosThetaStar",program_options::value<int>(),"")
          ("cosThetaStarCS",program_options::value<int>(),"")
          ("deltaEtaJets",program_options::value<int>(),"")
          ("productEtaJets",program_options::value<int>(),"")
          ("nJetsInRapidityGap",program_options::value<int>(),"")
          
          ("deltaPhiJets",program_options::value<int>(),"")
          ("deltaRJets",program_options::value<int>(),"")
          ("deltaEtaMuons",program_options::value<int>(),"")
          ("deltaPhiMuons",program_options::value<int>(),"")
          ("deltaRMuons",program_options::value<int>(),"")
      
          ("relIsoMu1",program_options::value<int>(),"")
          ("relIsoMu2",program_options::value<int>(),"")
          ("ht",program_options::value<int>(),"")
          ("nJets",program_options::value<int>(),"")
          ("htInRapidityGap",program_options::value<int>(),"")
          ("weightsDirName",program_options::value<std::string>(),"")
          ("vbfSelection",program_options::value<int>(),"")

          ("nVtx",program_options::value<int>(),"")
  
          ("sigFile",program_options::value<std::vector<std::string> >(),"")
          ("bakFile",program_options::value<std::vector<std::string> >(),"")
          ("sigWeight",program_options::value<std::vector<float> >(),"")
          ("bakWeight",program_options::value<std::vector<float> >(),"")

          ("nTrees",program_options::value<unsigned>(),"")
          ("nEventsMin",program_options::value<unsigned>(),"")
      ;
    
      program_options::variables_map optionMap;
      program_options::store(program_options::parse_config_file(infile, optionDesc), optionMap);
      program_options::notify(optionMap);    
      infile.close();
  
      TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");
    
      if (optionMap.count("mDiMu") && optionMap["mDiMu"].as<int>() == 1)
        reader->AddVariable("mDiMu",&mDiMu);
      else
        reader->AddSpectator("mDiMu",&mDiMu);
  
      if (optionMap.count("ptDiMu") && optionMap["ptDiMu"].as<int>() == 1)
        reader->AddVariable("ptDiMu",&ptDiMu);
      else
        reader->AddSpectator("ptDiMu",&ptDiMu);
  
      if (optionMap.count("yDiMu") && optionMap["yDiMu"].as<int>() == 1)
        reader->AddVariable("yDiMu",&yDiMu);
      else
        reader->AddSpectator("yDiMu",&yDiMu);
  
      if (optionMap.count("mDiJet") && optionMap["mDiJet"].as<int>() == 1)
        reader->AddVariable("mDiJet",&mDiJet);
      else
        reader->AddSpectator("mDiJet",&mDiJet);
  
      if (optionMap.count("ptDiJet") && optionMap["ptDiJet"].as<int>() == 1)
        reader->AddVariable("ptDiJet",&ptDiJet);
      else
        reader->AddSpectator("ptDiJet",&ptDiJet);
  
      if (optionMap.count("yDiJet") && optionMap["yDiJet"].as<int>() == 1)
        reader->AddVariable("yDiJet",&yDiJet);
      else
        reader->AddSpectator("yDiJet",&yDiJet);
  
      if (optionMap.count("ptMu1") && optionMap["ptMu1"].as<int>() == 1)
        reader->AddVariable("ptMu1",&ptMu1);
      else
        reader->AddSpectator("ptMu1",&ptMu1);
  
      if (optionMap.count("ptMu2") && optionMap["ptMu2"].as<int>() == 1)
        reader->AddVariable("ptMu2",&ptMu2);
      else
        reader->AddSpectator("ptMu2",&ptMu2);
  
      if (optionMap.count("etaMu1") && optionMap["etaMu1"].as<int>() == 1)
        reader->AddVariable("etaMu1",&etaMu1);
      else
        reader->AddSpectator("etaMu1",&etaMu1);
  
      if (optionMap.count("etaMu2") && optionMap["etaMu2"].as<int>() == 1)
        reader->AddVariable("etaMu2",&etaMu2);
      else
        reader->AddSpectator("etaMu2",&etaMu2);
  
      if (optionMap.count("ptJet1") && optionMap["ptJet1"].as<int>() == 1)
        reader->AddVariable("ptJet1",&ptJet1);
      else
        reader->AddSpectator("ptJet1",&ptJet1);
  
      if (optionMap.count("ptJet2") && optionMap["ptJet2"].as<int>() == 1)
        reader->AddVariable("ptJet2",&ptJet2);
      else
        reader->AddSpectator("ptJet2",&ptJet2);
  
      if (optionMap.count("etaJet1") && optionMap["etaJet1"].as<int>() == 1)
        reader->AddVariable("etaJet1",&etaJet1);
      else
        reader->AddSpectator("etaJet1",&etaJet1);
  
      if (optionMap.count("etaJet2") && optionMap["etaJet2"].as<int>() == 1)
        reader->AddVariable("etaJet2",&etaJet2);
      else
        reader->AddSpectator("etaJet2",&etaJet2);
  
      if (optionMap.count("cosThetaStar") && optionMap["cosThetaStar"].as<int>() == 1)
        reader->AddVariable("cosThetaStar",&cosThetaStar);
      else
        reader->AddSpectator("cosThetaStar",&cosThetaStar);
  
      if (optionMap.count("cosThetaStarCS") && optionMap["cosThetaStarCS"].as<int>() == 1)
        reader->AddVariable("cosThetaStarCS",&cosThetaStarCS);
      else
        reader->AddSpectator("cosThetaStarCS",&cosThetaStarCS);
  
      if (optionMap.count("deltaEtaJets") && optionMap["deltaEtaJets"].as<int>() == 1)
        reader->AddVariable("deltaEtaJets",&deltaEtaJets);
      else
        reader->AddSpectator("deltaEtaJets",&deltaEtaJets);
  
      if (optionMap.count("productEtaJets") && optionMap["productEtaJets"].as<int>() == 1)
        reader->AddVariable("productEtaJets",&productEtaJets);
      else
        reader->AddSpectator("productEtaJets",&productEtaJets);
  
      if (optionMap.count("nJetsInRapidityGap") && optionMap["nJetsInRapidityGap"].as<int>() == 1)
        reader->AddVariable("nJetsInRapidityGap",&nJetsInRapidityGap);
      else
        reader->AddSpectator("nJetsInRapidityGap",&nJetsInRapidityGap);
  
      if (optionMap.count("deltaEtaMuons") && optionMap["deltaEtaMuons"].as<int>() == 1)
        reader->AddVariable("deltaEtaMuons",&deltaEtaMuons);
      else
        reader->AddSpectator("deltaEtaMuons",&deltaEtaMuons);
  
      if (optionMap.count("deltaPhiMuons") && optionMap["deltaPhiMuons"].as<int>() == 1)
        reader->AddVariable("deltaPhiMuons",&deltaPhiMuons);
      else
        reader->AddSpectator("deltaPhiMuons",&deltaPhiMuons);
  
      if (optionMap.count("deltaRMuons") && optionMap["deltaRMuons"].as<int>() == 1)
        reader->AddVariable("deltaRMuons",&deltaRMuons);
      else
        reader->AddSpectator("deltaRMuons",&deltaRMuons);
  
      if (optionMap.count("deltaPhiJets") && optionMap["deltaPhiJets"].as<int>() == 1)
        reader->AddVariable("deltaPhiJets",&deltaPhiJets);
      else
        reader->AddSpectator("deltaPhiJets",&deltaPhiJets);
  
      if (optionMap.count("deltaRJets") && optionMap["deltaRJets"].as<int>() == 1)
        reader->AddVariable("deltaRJets",&deltaRJets);
      else
        reader->AddSpectator("deltaRJets",&deltaRJets);
  
      if (optionMap.count("relIsoMu1") && optionMap["relIsoMu1"].as<int>() == 1)
        reader->AddVariable("relIsoMu1",&relIsoMu1);
      else
        reader->AddSpectator("relIsoMu1",&relIsoMu1);
  
      if (optionMap.count("relIsoMu2") && optionMap["relIsoMu2"].as<int>() == 1)
        reader->AddVariable("relIsoMu2",&relIsoMu2);
      else
        reader->AddSpectator("relIsoMu2",&relIsoMu2);
  
      if (optionMap.count("ht") && optionMap["ht"].as<int>() == 1)
        reader->AddVariable("ht",&ht);
      else
        reader->AddSpectator("ht",&ht);
  
      if (optionMap.count("nJets") && optionMap["nJets"].as<int>() == 1)
        reader->AddVariable("nJets",&nJets);
      else
        reader->AddSpectator("nJets",&nJets);

      if (optionMap.count("nVtx") && optionMap["nVtx"].as<int>() == 1)
        reader->AddVariable("nVtx",&nVtx);
      else
        reader->AddSpectator("nVtx",&nVtx);
  
      if (optionMap.count("htInRapidityGap") && optionMap["htInRapidityGap"].as<int>() == 1)
        reader->AddVariable("htInRapidityGap",&htInRapidityGap);
      else
        reader->AddSpectator("htInRapidityGap",&htInRapidityGap);
    
      std::string weightsDirName;
      if (optionMap.count("weightsDirName"))
          weightsDirName = optionMap["weightsDirName"].as<std::string>();
      else
      {
        std::cout << "Error: Config File: " << *configFileName
              << " does not have option 'weightsDirName', it is required, exiting." << std::endl;
        throw;
      }
  
      reader->BookMVA("BDT",std::string(weightsDirName).append("/TMVAClassification_BDT.weights.xml").c_str());
      reader->BookMVA("Likelihood",std::string(weightsDirName).append("/TMVAClassification_Likelihood.weights.xml").c_str());
      readers_.insert(make_pair(*configFileName,reader));
    } // end Config file loop
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
  std::map<std::string,TMVA::Reader*>::const_iterator reader = readers_.begin();
  std::map<std::string,TMVA::Reader*>::const_iterator readerEnd = readers_.end();
  for(;reader != readerEnd; reader++)
  {
    delete reader->second;;
  }
}

float
MVA::getMVA(const std::string configFileName, const std::string mvaName)
{
  float value = -1000.0;
  if(readers_.count(configFileName))
  {
    value = readers_[configFileName]->EvaluateMVA(mvaName.c_str());
  }
  return value;
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
  cosThetaStarCS=-10.0;
  deltaEtaJets=-10.0;
  productEtaJets=0.0;
  nJetsInRapidityGap=0.0;

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

  nVtx=0.0;

  weight=1.0;
}
