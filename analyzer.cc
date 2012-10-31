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
#include <TRandom3.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "DataFormats.h"
#include "helpers.h"
#include "mva.h"
#include "LumiReweightingStandAlone.h"

#include "boost/program_options.hpp"
#include "boost/regex.hpp"

#include <boost/lexical_cast.hpp>
//Defines method of std::string that appends any type :-)
#define appendAny(a) append(boost::lexical_cast<std::string>(a))

#include <limits.h>

#define JETPUID
#define PUREWEIGHT
//#define SMEARING

using namespace std;
using namespace boost;

void fillMuonHist(TH1F* hist, _MuonInfo& mu1, _MuonInfo& mu2);
void printStationMiss(_MuonInfo& mu1, _MuonInfo& mu2, _EventInfo& eventInfo, std::string & testString, unsigned & testCounter);

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  const char* optionIntro = "H->MuMu Analyzer\n\nUsage: ./analyzer [--help] [--train] [--maxEvents N] <outputFileName.root> <inputFileName.root> [<inputFileName2.root>...]\n\nAllowed Options";
  program_options::options_description optionDesc(optionIntro);
  optionDesc.add_options()
      ("help,h", "Produce Help Message")
      ("trainingTree,t",program_options::value<string>(), "Create Training Tree File with filename")
      ("filenames",program_options::value<vector<string> >(), "Input & Output File Names, put output name first followed by all input file names")
      ("maxEvents,m",program_options::value<int>(), "Maximum Number of Events to Process")
  ;
  
  program_options::positional_options_description optionPos;
  optionPos.add("filenames",-1);
  
  program_options::variables_map optionMap;
  program_options::store(program_options::command_line_parser(argc, argv).options(optionDesc).positional(optionPos).run(), optionMap);
  program_options::notify(optionMap);    
  
  if (optionMap.count("help")) 
  {
      cout << optionDesc << "\n";
      return 1;
  }

  std::vector<std::string> filenames;
  vector<string>::const_iterator filename;
  std::string outputFileName;
  if (optionMap.count("filenames")>0)
  {
     filenames = optionMap["filenames"].as<vector<string> >();
     if(filenames.size()<2)
     {
       cout << "Error: Need both input file and output file names, exiting." << endl;
       return 1;
     }
     outputFileName = filenames[0];
     filenames.erase(filenames.begin());
  }
  else
  {
     cout << "Error: Input file name  and ouput file name arguments required, exiting." << endl;
     return 1;
  }

  int maxEvents = std::numeric_limits<int>::max();
  if (optionMap.count("maxEvents")) 
  {
      int tmp = optionMap["maxEvents"].as<int>();
      if (tmp > 0)
      {
        maxEvents = tmp;
      }
  }
  cout << "maxEvents = "<< maxEvents << "\n";

  /////////////////////////////
  //////////// Setup //////////
  /////////////////////////////

  float minMmm = 70.0;
  float maxMmm = 200.0;
  //float minMmm = 110.0;
  //float maxMmm = 150.0;

  float minBlind = 120;
  float maxBlind = 130;

  float calib = -0.1;
  float calibSysSmear = 0.2;
  float resSmear = 1.169; // should all be around 1; a ratio of resolutions
  float resSysSmear = 0.2; // Error on that ratio

  std::vector<int> allowedHLTPaths;
  allowedHLTPaths.push_back(0); //IsoMu24

  // Check to see if it is data
  bool isData = false;
  std::vector<std::string> dataWords;
  dataWords.push_back("Run2012");
  dataWords.push_back("Run2011");
  dataWords.push_back("SingleMu");
  dataWords.push_back("DoubleMu");
  std::vector<std::string>::const_iterator dataWord;
  for(dataWord = dataWords.begin(); dataWord != dataWords.end();dataWord++)
  {
    regex re(*dataWord);
    if(regex_search(filenames[0],re))
    {
        isData = true;
    }
  }

  if(isData)
    std::cout << "This is a Real Data Sample\n";
  else
    std::cout << "This is a MC Sample\n";

  ////////////
  
  TChain * tree = new TChain("tree");

  cout << "Input File Names: \n"; 
  for(filename = filenames.begin();filename != filenames.end();filename++)
  {
    cout<<"  "<< *filename << endl;
    tree->AddFile(filename->c_str());
  }

  cout << "Output File Name: " << outputFileName << endl;
  TFile * outFile = new TFile(outputFileName.c_str(),"RECREATE");

  std::string trainingTreeFileName = "";
  if (optionMap.count("trainingTree")) 
  {
      cout << "Training enabled" << "\n";
      trainingTreeFileName = optionMap["trainingTree"].as<string>();
      cout << "Training Tree File Name: " << trainingTreeFileName << "\n";
  }

  //////////////////////////
  // Tree Branches

  _MuonInfo reco1, reco2;

  tree->SetBranchAddress("reco1", &reco1);
  tree->SetBranchAddress("reco2", &reco2);

  float recoCandMass, recoCandPt, recoCandY, recoCandPhi;

  tree->SetBranchAddress("recoCandMass", &recoCandMass);
  tree->SetBranchAddress("recoCandPt"  , &recoCandPt );
  tree->SetBranchAddress("recoCandY"  , &recoCandY );
  tree->SetBranchAddress("recoCandPhi"  , &recoCandPhi );

  float trueMass=-99999.0;
  if(!isData && tree->GetBranchStatus("trueMass"))
    tree->SetBranchAddress("trueMass", &trueMass);

  _PFJetInfo jets;
  tree->SetBranchAddress("pfJets",&jets);

#ifdef JETPUID
  float puJetFullDisc[10];
  float puJetSimpleDisc[10];
  float puJetCutDisc[10];

  tree->SetBranchAddress("puJetFullDisc",&puJetFullDisc);
  tree->SetBranchAddress("puJetSimpleDisc",&puJetSimpleDisc);
  tree->SetBranchAddress("puJetCutDisc",&puJetCutDisc);

  int puJetFullId[10];
  int puJetSimpleId[10];
  int puJetCutId[10];

  tree->SetBranchAddress("puJetFullId",&puJetFullId);
  tree->SetBranchAddress("puJetSimpleId",&puJetSimpleId);
  tree->SetBranchAddress("puJetCutId",&puJetCutId);
#endif

  int nPU=0;
#ifdef PUREWEIGHT
  if (!isData)
  {
    tree->SetBranchAddress("nPU",&nPU);
  }
#endif
  _VertexInfo vertexInfo;
  tree->SetBranchAddress("vertexInfo",&vertexInfo);
  _EventInfo eventInfo;
  tree->SetBranchAddress("eventInfo",&eventInfo);
  _MetInfo met;
  tree->SetBranchAddress("met",&met);

  //////////////////////////
  // Histograms
  std::map<std::string,TH1F*> histMap;
  std::map<std::string,TH2F*> histMap2D;
  std::map<std::string,TH1F*>::iterator histMapIter;
  std::map<std::string,TH2F*>::iterator histMap2DIter;

  TH1F* mDiMu = new TH1F("mDiMu","DiMuon Mass",1600,0,400);
  histMap.insert(make_pair("mDiMu",mDiMu));

  TH1F* mDiJet = new TH1F("mDiJet","DiJet Mass",500,0,2000);
  histMap.insert(make_pair("mDiJet",mDiJet));

  TH1F* ptDiMu = new TH1F("ptDiMu","DiMuon Pt",250,0,500);
  histMap.insert(make_pair("ptDiMu",ptDiMu));

  TH1F* yDiMu = new TH1F("yDiMu","DiMuon Rapidity",100,-4,4);
  histMap.insert(make_pair("yDiMu",yDiMu));

  TH2F* yVptDiMu = new TH2F("yVptDiMu","DiMuon Rapidity v. p_{T}",250,0,500,100,0,4);
  histMap2D.insert(make_pair("yVptDiMu",yVptDiMu));
  TH2F* ptVmDiMu = new TH2F("ptVmDiMu","DiMuon p_{T} v. Mass",1600,0,400,250,0,250);
  histMap2D.insert(make_pair("ptVmDiMu",ptVmDiMu));
  TH2F* yVmDiMu = new TH2F("yVmDiMu","DiMuon |y| v. Mass",1600,0,400,100,0,4);
  TH2F* phiVmDiMu = new TH2F("phiVmDiMu","DiMuon #phi v. Mass",1600,0,400,100,0,3.2);
  histMap2D.insert(make_pair("phiVmDiMu",phiVmDiMu));

  TH1F* ptMu1 = new TH1F("ptMu1","Leading Muon Pt",100,0,200);
  histMap.insert(make_pair("ptMu1",ptMu1));
  TH1F* ptMu2 = new TH1F("ptMu2","Sub-Leading Muon Pt",100,0,200);
  histMap.insert(make_pair("ptMu2",ptMu2));
  TH1F* ptJet1 = new TH1F("ptJet1","Leading Jet Pt",200,0,1000);
  histMap.insert(make_pair("ptJet1",ptJet1));
  TH1F* ptJet2 = new TH1F("ptJet2","Sub-Leading Jet Pt",200,0,1000);
  histMap.insert(make_pair("ptJet2",ptJet2));

  TH1F* etaMu1 = new TH1F("etaMu1","Leading Muon #eta",50,-2.5,2.5);
  histMap.insert(make_pair("etaMu1",etaMu1));
  TH1F* etaMu2 = new TH1F("etaMu2","Sub-Leading Muon #eta",50,-2.5,2.5);
  histMap.insert(make_pair("etaMu2",etaMu2));
  TH1F* etaJet1 = new TH1F("etaJet1","Leading Jet #eta",50,-5.0,5.0);
  histMap.insert(make_pair("etaJet1",etaJet1));
  TH1F* etaJet2 = new TH1F("etaJet2","Sub-Leading Jet #eta",50,-5.0,5.0);
  histMap.insert(make_pair("etaJet2",etaJet2));

  TH1F* deltaEtaJetsHist = new TH1F("deltaEtaJets","#Delta#eta Jets",50,0.0,10.0);
  histMap.insert(make_pair("deltaEtaJets",deltaEtaJetsHist));
  TH1F* deltaPhiJetsHist = new TH1F("deltaPhiJets","#Delta#phi Jets",50,0.0,3.2);
  histMap.insert(make_pair("deltaPhiJets",deltaPhiJetsHist));
  TH1F* deltaRJetsHist = new TH1F("deltaRJets","#Delta R Jets",50,0.0,10.0);
  histMap.insert(make_pair("deltaRJets",deltaRJetsHist));

  TH1F* deltaEtaMuonsHist = new TH1F("deltaEtaMuons","#Delta#eta Jets",50,0.0,10.0);
  histMap.insert(make_pair("deltaEtaMuons",deltaEtaMuonsHist));
  TH1F* deltaPhiMuonsHist = new TH1F("deltaPhiMuons","#Delta#phi Jets",50,0.0,3.2);
  histMap.insert(make_pair("deltaPhiMuons",deltaPhiMuonsHist));
  TH1F* deltaRMuonsHist = new TH1F("deltaRMuons","#Delta R Jets",50,0.0,10.0);
  histMap.insert(make_pair("deltaRMuons",deltaRMuonsHist));

  TH1F* countsHist = new TH1F("countsHist","Event Counts",10,0.0,10.0);
  countsHist->GetXaxis()->SetBinLabel(1,"total");
  countsHist->GetXaxis()->SetBinLabel(2,"2#mu ID");
  countsHist->GetXaxis()->SetBinLabel(3,"HLT");
  countsHist->GetXaxis()->SetBinLabel(4,"Charge");
  countsHist->GetXaxis()->SetBinLabel(5,"m_{#mu#mu}");
  countsHist->GetXaxis()->SetBinLabel(6,"Inc Pre");
  countsHist->GetXaxis()->SetBinLabel(7,"VBF Pre");
  histMap.insert(make_pair("countsHist",countsHist));

  TH1F* countsHist2 = new TH1F("countsHist2","Event Counts",14,0.0,14.0);
  countsHist2->GetXaxis()->SetBinLabel(1,"total");
  countsHist2->GetXaxis()->SetBinLabel(2,"total");
  countsHist2->GetXaxis()->SetBinLabel(3,"global");
  countsHist2->GetXaxis()->SetBinLabel(4,"PF");
  countsHist2->GetXaxis()->SetBinLabel(5,"pt");
  countsHist2->GetXaxis()->SetBinLabel(6,"eta");
  countsHist2->GetXaxis()->SetBinLabel(7,"tracker");
  countsHist2->GetXaxis()->SetBinLabel(8,"iso");
  countsHist2->GetXaxis()->SetBinLabel(9,"d0");
  countsHist2->GetXaxis()->SetBinLabel(10,"dz");
  countsHist2->GetXaxis()->SetBinLabel(11,"#mu hits");
  countsHist2->GetXaxis()->SetBinLabel(12,"pixel");
  countsHist2->GetXaxis()->SetBinLabel(13,"stations");
  countsHist2->GetXaxis()->SetBinLabel(14,"#chi^2");
  histMap.insert(make_pair("countsHist2",countsHist2));

  TH1F* cosThetaStarHist = new TH1F("cosThetaStar","cos(#theta^{*})",50,-1.,1.);
  histMap.insert(make_pair("cosThetaStar",cosThetaStarHist));
  TH1F* cosThetaStarCSHist = new TH1F("cosThetaStarCS","cos(#theta^{*}_{CS})",50,-1.,1.);
  histMap.insert(make_pair("cosThetaStarCS",cosThetaStarCSHist));

  TH1F* puJetIDSimpleDiscJet1Hist = new TH1F("puJetIDSimpleDiscJet1","PU Jet ID--Simple Discriminator Leading Jet",50,-1.,1.);
  histMap.insert(make_pair("puJetIDSimpleDiscJet1",puJetIDSimpleDiscJet1Hist));
  TH1F* puJetIDSimpleDiscJet2Hist = new TH1F("puJetIDSimpleDiscJet2","PU Jet ID--Simple Discriminator Sub-Leading Jet",50,-1.,1.);
  histMap.insert(make_pair("puJetIDSimpleDiscJet2",puJetIDSimpleDiscJet2Hist));
  TH1F* puJetIDSimpleDiscJet3Hist = new TH1F("puJetIDSimpleDiscJet3","PU Jet ID--Simple Discriminator 3rd Leading Jet",50,-1.,1.);
  histMap.insert(make_pair("puJetIDSimpleDiscJet3",puJetIDSimpleDiscJet3Hist));

  TH1F* puJetIDSimpleJet1Hist = new TH1F("puJetIDSimpleJet1","PU Jet ID--Simple Loose Leading Jet",2,0,2);
  histMap.insert(make_pair("puJetIDSimpleJet1",puJetIDSimpleJet1Hist));
  puJetIDSimpleJet1Hist->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet1Hist->GetXaxis()->SetBinLabel(2,"Pass");
  TH1F* puJetIDSimpleJet2Hist = new TH1F("puJetIDSimpleJet2","PU Jet ID--Simple Loose Sub-Leading Jet",2,-0,2);
  histMap.insert(make_pair("puJetIDSimpleJet2",puJetIDSimpleJet2Hist));
  puJetIDSimpleJet2Hist->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet2Hist->GetXaxis()->SetBinLabel(2,"Pass");
  TH1F* puJetIDSimpleJet3Hist = new TH1F("puJetIDSimpleJet3","PU Jet ID--Simple Loose 3rd Leading Jet",2,-0,2);
  histMap.insert(make_pair("puJetIDSimpleJet3",puJetIDSimpleJet3Hist));
  puJetIDSimpleJet3Hist->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet3Hist->GetXaxis()->SetBinLabel(2,"Pass");

  TH1F* BDTHistMuonOnly = new TH1F("BDTHistMuonOnly","BDT Discriminator",2000,-1,1);
  histMap.insert(make_pair("BDTHistMuonOnly",BDTHistMuonOnly));
  TH1F* likelihoodHistMuonOnly = new TH1F("likelihoodHistMuonOnly","Likelihood Discriminator",2000,-1,1);
  histMap.insert(make_pair("likelihoodHistMuonOnly",likelihoodHistMuonOnly));

  TH1F* BDTHistVBF = new TH1F("BDTHistVBF","BDT Discriminator",2000,-1,1);
  histMap.insert(make_pair("BDTHistVBF",BDTHistVBF));
  TH1F* likelihoodHistVBF = new TH1F("likelihoodHistVBF","Likelihood Discriminator",2000,-1,1);
  histMap.insert(make_pair("likelihoodHistVBF",likelihoodHistVBF));

  TH2F* BDTHistMuonOnlyVMass = new TH2F("BDTHistMuonOnlyVMass","BDT Discriminator",1600,0,400,2000,-1,1);
  histMap2D.insert(make_pair("BDTHistMuonOnlyVMass",BDTHistMuonOnlyVMass));
  TH2F* likelihoodHistMuonOnlyVMass = new TH2F("likelihoodHistMuonOnlyVMass","Likelihood Discriminator",1600,0,400,2000,-1,1);
  histMap2D.insert(make_pair("likelihoodHistMuonOnlyVMass",likelihoodHistMuonOnlyVMass));
  TH2F* BDTHistVBFVMass = new TH2F("BDTHistVBFVMass","BDT Discriminator",1600,0,400,2000,-1,1);
  histMap2D.insert(make_pair("BDTHistVBFVMass",BDTHistVBFVMass));
  TH2F* likelihoodHistVBFVMass = new TH2F("likelihoodHistVBFVMass","Likelihood Discriminator",1600,0,400,2000,-1,1);
  histMap2D.insert(make_pair("likelihoodHistVBFVMass",likelihoodHistVBFVMass));

  TH1F* relIsoMu1Hist = new TH1F("relIsoMu1","",1000,0,10.0);
  histMap.insert(make_pair("relIsoMu1",relIsoMu1Hist));
  TH1F* relIsoMu2Hist = new TH1F("relIsoMu2","",1000,0,10.0);
  histMap.insert(make_pair("relIsoMu2",relIsoMu2Hist));

  TH1F* nJetsHist = new TH1F("nJets","",11,0,11);
  histMap.insert(make_pair("nJets",nJetsHist));
  TH1F* htHist = new TH1F("ht","",200,0,2000);
  histMap.insert(make_pair("ht",htHist));
  TH1F* nJetsInRapidityGapHist = new TH1F("nJetsInRapidityGap","",11,0,11);
  histMap.insert(make_pair("nJetsInRapidityGap",nJetsInRapidityGapHist));
  TH1F* htInRapidityGapHist = new TH1F("htInRapidityGap","",200,0,2000);
  histMap.insert(make_pair("htInRapidityGap",htInRapidityGapHist));

  TH1F* nPUHist = new TH1F("nPU","",100,0,100);
  histMap.insert(make_pair("nPU",nPUHist));
  TH1F* nVtxHist = new TH1F("nVtx","",100,0,100);
  histMap.insert(make_pair("nVtx",nVtxHist));
  TH1F* metHist = new TH1F("met","",160,0,800);
  histMap.insert(make_pair("met",metHist));
  TH1F* weightHist = new TH1F("weight","",500,0,5.0);
  histMap.insert(make_pair("weight",weightHist));

  for(histMapIter = histMap.begin(); histMapIter != histMap.end(); histMapIter++)
  {
    histMapIter->second->Sumw2();
  }
  for(histMap2DIter = histMap2D.begin(); histMap2DIter != histMap2D.end(); histMap2DIter++)
  {
    histMap2DIter->second->Sumw2();
  }

  // Other Sets of Hists
  std::map<std::string,TH1F*> histMap4GeVWindow;
  std::map<std::string,TH1F*> histMapPtDiMu100;
  std::map<std::string,TH1F*> histMapVBFPresel;
  std::map<std::string,TH1F*> histMapIncPresel;
  std::map<std::string,TH1F*> histMapNotBlindWindow;
  std::map<std::string,TH1F*> histMapUpperControlRegion;
  std::map<std::string,TH1F*> histMapLowerControlRegion;
  std::map<std::string,TH1F*> histMapVBFLoose;
  std::map<std::string,TH1F*> histMapVBFMedium;
  std::map<std::string,TH1F*> histMapVBFTight;
  std::map<std::string,TH1F*> histMapVBFVeryTight;
  std::map<std::string,TH1F*> histMapPt0to30;
  std::map<std::string,TH1F*> histMapPt30to50;
  std::map<std::string,TH1F*> histMapPt50to125;
  std::map<std::string,TH1F*> histMapPt125to250;
  std::map<std::string,TH1F*> histMapPt250;
  std::map<std::string,TH1F*> histMapIncBDTSig80;
  std::map<std::string,TH1F*> histMapVBFBDTSig80;
  for(histMapIter = histMap.begin(); histMapIter != histMap.end(); histMapIter++)
  {
    TH1F* tmp;
    tmp = (TH1F*) histMapIter->second->Clone();
    histMap4GeVWindow.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapPtDiMu100.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapVBFPresel.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapIncPresel.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapNotBlindWindow.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapUpperControlRegion.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapLowerControlRegion.insert(make_pair(histMapIter->first,tmp));

    tmp = (TH1F*) histMapIter->second->Clone();
    histMapVBFLoose.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapVBFMedium.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapVBFTight.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapVBFVeryTight.insert(make_pair(histMapIter->first,tmp));

    tmp = (TH1F*) histMapIter->second->Clone();
    histMapPt0to30.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapPt30to50.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapPt50to125.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapPt125to250.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapPt250.insert(make_pair(histMapIter->first,tmp));

    tmp = (TH1F*) histMapIter->second->Clone();
    histMapIncBDTSig80.insert(make_pair(histMapIter->first,tmp));
    tmp = (TH1F*) histMapIter->second->Clone();
    histMapVBFBDTSig80.insert(make_pair(histMapIter->first,tmp));
  }

  std::map<std::string,TH2F*> histMap2DCalibUp;
  std::map<std::string,TH2F*> histMap2DCalibDown;
  std::map<std::string,TH2F*> histMap2DResUp;
  std::map<std::string,TH2F*> histMap2DResDown;
  for(histMap2DIter = histMap2D.begin(); histMap2DIter != histMap2D.end(); histMap2DIter++)
  {
    TH2F* tmp;
    tmp = (TH2F*) histMap2DIter->second->Clone();
    histMap2DCalibUp.insert(make_pair(histMap2DIter->first,tmp));
    tmp = (TH2F*) histMap2DIter->second->Clone();
    histMap2DCalibDown.insert(make_pair(histMap2DIter->first,tmp));

    tmp = (TH2F*) histMap2DIter->second->Clone();
    histMap2DResUp.insert(make_pair(histMap2DIter->first,tmp));
    tmp = (TH2F*) histMap2DIter->second->Clone();
    histMap2DResDown.insert(make_pair(histMap2DIter->first,tmp));
  }

  //////////////////////////
  //for MVA

  std::vector<std::string> mvaConfigNames;
  mvaConfigNames.push_back("inclusive.cfg");
  mvaConfigNames.push_back("vbf.cfg");
  MVA mva(mvaConfigNames,trainingTreeFileName);

  TRandom3 random(1457);

  //////////////////////////
  //for PU reweighting

#ifdef PUREWEIGHT
  reweight::LumiReWeighting lumiWeights("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root","pileupDists/PileUpHist2012A.root","pileup","pileup");
#endif

  const double SQRT2 = sqrt(2);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  unsigned nEvents = tree->GetEntries();
  cout << "nEvents: " << nEvents << endl;

  unsigned reportEach=1000;
  if (nEvents/1000>reportEach)
    reportEach = nEvents/1000;

  unsigned testCounter = 0;
  string testString;
  
  for(unsigned i=0; i<nEvents;i++)
  {
    if(i >= maxEvents)
      continue;
    tree->GetEvent(i);
    if (i % reportEach == 0) cout << "Event: " << i << endl;

    fillMuonHist(countsHist2, reco1, reco2);
    //printStationMiss(reco1,reco2,eventInfo,testString,testCounter);

    mva.resetValues();
    mva.mDiMu = recoCandMass;
#ifdef SMEARING
    if(!isData)
    {
      mva.mDiMu = smearMC(trueMass,recoCandMass,calib,resSmear,random);
    }
#endif
    bool inBlindWindow = mva.mDiMu < maxBlind && mva.mDiMu > minBlind;

    double weight = 1.0;
#ifdef PUREWEIGHT
    if (!isData)
    {
      weight = lumiWeights.weight(nPU);
    }
#endif

    countsHist->Fill(0.0, weight);

    if (!isKinTight_2012(reco1) || !isKinTight_2012(reco2))
        continue;

    countsHist->Fill(1.0, weight);

    if (!isHltMatched(reco1,reco2,allowedHLTPaths) && isData)
        continue;

    countsHist->Fill(2.0, weight);

    if (reco1.charge*reco2.charge != -1)
        continue;

    countsHist->Fill(3.0, weight);

    if (mva.mDiMu < minMmm || mva.mDiMu > maxMmm)
        continue;

    countsHist->Fill(4.0, weight);

    _MuonInfo muon1;
    _MuonInfo muon2;
    if(reco1.pt>reco2.pt)
    {
        muon1 = reco1;
        muon2 = reco2;
    }
    else
    {
        muon1 = reco2;
        muon2 = reco1;
    }

    mva.weight = weight;

    mva.ptMu1=muon1.pt;
    mva.ptMu2=muon2.pt;
    mva.etaMu1=muon1.eta;
    mva.etaMu2=muon2.eta;
    mva.deltaEtaMuons=fabs(muon1.eta-muon2.eta);
    mva.relIsoMu1 = getRelIso(muon1);
    mva.relIsoMu2 = getRelIso(muon2);

    mva.ptDiMu = recoCandPt;
    mva.yDiMu = recoCandY;

    float mDiMuCalibUp = mva.mDiMu+calibSysSmear;
    float mDiMuCalibDown = mva.mDiMu-calibSysSmear;

    float mDiMuResUp = smearMC(trueMass,recoCandMass,calib,resSmear+resSysSmear,random);
    float mDiMuResDown = smearMC(trueMass,recoCandMass,calib,resSmear-resSysSmear,random);

    //////////////////////////////////////////
    //Computing CosTheta*

    TLorentzVector pMuon1;
    TLorentzVector pMuon2;
    pMuon1.SetPtEtaPhiM(muon1.pt,muon1.eta,muon1.phi,0.105);
    pMuon2.SetPtEtaPhiM(muon2.pt,muon2.eta,muon2.phi,0.105);
    TLorentzVector diMuon = pMuon1+pMuon2;

    mva.deltaPhiMuons = pMuon1.DeltaPhi(pMuon2);
    mva.deltaRMuons = pMuon1.DeltaR(pMuon2);

    TLorentzVector starMuon1 = pMuon1;
    TLorentzVector starMuon2 = pMuon2;
    TVector3 boost = diMuon.BoostVector();
    starMuon1.Boost(-boost);
    starMuon2.Boost(-boost);

    //if (muon1.charge>0)
    //std::cout << "Run: " << eventInfo.run << " lumi: " << eventInfo.lumi << " event: " << eventInfo.event << std::endl;
    if ((int) (eventInfo.event) % 2 == 0)
    {
        TVector3 directionOfBoost = starMuon1.BoostVector();
        mva.cosThetaStar = directionOfBoost.Dot(diMuon.BoostVector()) / (directionOfBoost.Mag()*diMuon.BoostVector().Mag());
    }
    else
    {
        TVector3 directionOfBoost = starMuon2.BoostVector();
        mva.cosThetaStar = directionOfBoost.Dot(diMuon.BoostVector()) / (directionOfBoost.Mag()*diMuon.BoostVector().Mag());
    }

    //////////////////////////////////////////
    //Computing CosTheta* Collins-Soper

    //std::cout << "muon1 charge: " << muon1.charge << "muon2 charge: "<<muon2.charge << std::endl;
    if (muon1.charge != muon2.charge)
    {
      // p1 is lepton
      // p2 is anti-lepton
      float p1Plus=-1e15;
      float p2Plus=-1e15;
      float p1Minus=-1e15;
      float p2Minus=-1e15;
      if (muon1.charge < 0)
      {
        p1Plus  = (pMuon1.E()+pMuon1.Pz())/SQRT2;
        p1Minus = (pMuon1.E()-pMuon1.Pz())/SQRT2;
        p2Plus  = (pMuon2.E()+pMuon2.Pz())/SQRT2;
        p2Minus = (pMuon2.E()-pMuon2.Pz())/SQRT2;
      }
      else
      {
        p1Plus  = (pMuon2.E()+pMuon2.Pz())/SQRT2;
        p1Minus = (pMuon2.E()-pMuon2.Pz())/SQRT2;
        p2Plus  = (pMuon1.E()+pMuon1.Pz())/SQRT2;
        p2Minus = (pMuon1.E()-pMuon1.Pz())/SQRT2;
      }
      mva.cosThetaStarCS = diMuon.Pz()/fabs(diMuon.Pz()) * 
                               2*(p1Plus*p2Minus-p1Minus*p2Plus) / 
                               (diMuon.Mag()*sqrt(diMuon.Mag2()+diMuon.Pt()*diMuon.Pt()));
    }

    // Computing nVtx Valid
    //for(unsigned iVtx=0;iVtx<vertexInfo.nVertices;iVtx++)
    //{
    //  if(vertexInfo.isValid[iVtx])
    //  {
    //    mva.nVtx++;
    //  }
    //}
    mva.nVtx = vertexInfo.nVertices;

    //////////////////////////////////////////
    // Filling Hists

#ifdef BLIND
    if (!(inBlindWindow && isData))
    {
#endif
    mDiMu->Fill(mva.mDiMu, weight);
    yVmDiMu->Fill(mva.mDiMu,fabs(mva.yDiMu), weight);
    ptVmDiMu->Fill(mva.mDiMu,mva.ptDiMu, weight);
    phiVmDiMu->Fill(mva.mDiMu,recoCandPhi, weight);
#ifdef BLIND
    }
#endif

    yDiMu->Fill(mva.yDiMu, weight);
    ptDiMu->Fill(mva.ptDiMu, weight);
    yVptDiMu->Fill(mva.ptDiMu,fabs(mva.yDiMu), weight);
    ptMu1->Fill(mva.ptMu1, weight);
    ptMu2->Fill(mva.ptMu2, weight);
    etaMu1->Fill(mva.etaMu1, weight);
    etaMu2->Fill(mva.etaMu2, weight);
    cosThetaStarHist->Fill(mva.cosThetaStar, weight);
    cosThetaStarCSHist->Fill(mva.cosThetaStarCS, weight);

    deltaPhiMuonsHist->Fill(mva.deltaPhiMuons, weight);
    deltaEtaMuonsHist->Fill(mva.deltaEtaMuons, weight);
    deltaRMuonsHist->Fill(mva.deltaRMuons, weight);

    relIsoMu1Hist->Fill(mva.relIsoMu1, weight);
    relIsoMu2Hist->Fill(mva.relIsoMu2, weight);

    nPUHist->Fill(nPU, weight);
    nVtxHist->Fill(mva.nVtx, weight);
    metHist->Fill(met.pt, weight);
    weightHist->Fill(weight);

    // Jet Part
    for(unsigned iJet=0; (iJet < jets.nJets && iJet < 10);iJet++)
    {
        if(jets.pt[iJet] > 30.0)
        {
          mva.nJets++;
          mva.ht += jets.pt[iJet];
        }
    }

    bool goodJets = false;
    if(jets.nJets>=2 && jets.pt[0]>30.0 && jets.pt[1]>30.0)
        goodJets = true;

//    if (mva.mDiMu > 140. && mva.mDiMu < 150. && goodJets)
//    {
//        testCounter++;
//        //std::cout <<eventInfo.run <<":"<<eventInfo.event <<"\n"<< std::endl;
//        testString.appendAny(eventInfo.run);
//        testString.append(":");
//        testString.appendAny(eventInfo.event);
//        testString.append("\n");
//    }


    if(goodJets)
    {
      TLorentzVector pJet1;
      TLorentzVector pJet2;
      pJet1.SetXYZM(jets.px[0],jets.py[0],jets.pz[0],jets.mass[0]);
      pJet2.SetXYZM(jets.px[1],jets.py[1],jets.pz[1],jets.mass[1]);
      TLorentzVector diJet = pJet1+pJet2;

      double dEtaJets = fabs(jets.eta[0]-jets.eta[1]);
      double etaJetProduct = jets.eta[0]*jets.eta[1];
      mva.deltaPhiJets = pJet1.DeltaPhi(pJet2);
      mva.deltaRJets = pJet1.DeltaR(pJet2);

      // Seeing if there are jets in the rapidity gap
      float etaMax = jets.eta[0];
      float etaMin = 9999999.0;
      if(etaMax < jets.eta[1])
      {
          etaMax = jets.eta[1];
          etaMin = jets.eta[0];
      }
      else
      {
          etaMin = jets.eta[1];
      }
      bool jetInRapidityGap=false;
      for(unsigned iJet=2; (iJet < jets.nJets && iJet < 10);iJet++)
      {
        if(jets.pt[iJet] > 30.0)
        {
          if(jets.eta[iJet] < etaMax && jets.eta[iJet] > etaMin)
          {
            jetInRapidityGap = true;
            mva.nJetsInRapidityGap++;
            mva.htInRapidityGap += jets.pt[iJet];
          }
        }
      }

#ifdef JETPUID
      for(unsigned iJet=0; (iJet < jets.nJets && iJet < 10);iJet++)
      {
        if(jets.pt[iJet]>30.0)
        {
          if (iJet==0)
            puJetIDSimpleDiscJet1Hist->Fill(puJetSimpleDisc[iJet], weight);
          else if (iJet==1)
            puJetIDSimpleDiscJet2Hist->Fill(puJetSimpleDisc[iJet], weight);
          else if (iJet==2)
            puJetIDSimpleDiscJet3Hist->Fill(puJetSimpleDisc[iJet], weight);

          if (iJet==0)
            puJetIDSimpleJet1Hist->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose), weight);
          else if (iJet==1)
            puJetIDSimpleJet2Hist->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose), weight);
          else if (iJet==2)
            puJetIDSimpleJet3Hist->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose), weight);
        }
      }
#endif

      mva.mDiJet = diJet.M();
      mva.yDiJet = diJet.Rapidity();
      mva.ptDiJet = diJet.Pt();
      mva.ptJet1 = pJet1.Pt();
      mva.ptJet2 = pJet2.Pt();
      mva.etaJet1 = pJet1.Eta();
      mva.etaJet2 = pJet2.Eta();
      mva.productEtaJets = etaJetProduct;
      mva.deltaEtaJets = dEtaJets;

      mDiJet->Fill(mva.mDiJet, weight);
      ptJet1->Fill(mva.ptJet1, weight);
      ptJet2->Fill(mva.ptJet2, weight);
      etaJet1->Fill(mva.etaJet1, weight);
      etaJet2->Fill(mva.etaJet2, weight);
      deltaEtaJetsHist->Fill(mva.deltaEtaJets, weight);
      deltaPhiJetsHist->Fill(mva.deltaPhiJets, weight);
      deltaRJetsHist->Fill(mva.deltaRJets, weight);
      nJetsInRapidityGapHist->Fill(mva.nJetsInRapidityGap, weight);
      htInRapidityGapHist->Fill(mva.htInRapidityGap, weight);
      nJetsHist->Fill(mva.nJets, weight);
      htHist->Fill(mva.ht, weight);
    }
  
//HIG-12-007 PAS H->tautau
//The VBF category requires at least two jets with pT > 30 GeV/c, |η1 − η2 | > 4.0 and
//η1 · η2 < 0 (with η1 the pseudorapidty of the leading jet and η2 the pseudorapidity
//of the subleading jet), and a di-jet invariant mass m12 > 400 GeV/c2 , with no other
//jet with pT > 30 GeV/c in the rapidity region between the two jets.

    mva.writeEvent();

    bool vbfPreselection = mva.mDiJet>300.0 && mva.deltaEtaJets>3.0 && mva.productEtaJets<0.0 && mva.nJetsInRapidityGap == 0;
    //if(vbfPreselection)
    //  std::cout << "VBF Preselected!!";

    if(vbfPreselection)
        countsHist->Fill(6);
    else
        countsHist->Fill(5);
    

    bool vbfVeryTight = false;
    bool vbfTight = false;
    bool vbfMedium = false;
    bool vbfLoose = false;
    if(vbfPreselection && mva.mDiJet>700.0 && mva.deltaEtaJets>5.)
    {
        vbfVeryTight=true;
    }
    else if(vbfPreselection && mva.mDiJet>400.0 && mva.deltaEtaJets>5.)
    {
        vbfTight=true;
    }
    else if(vbfPreselection && mva.mDiJet>400.0 && mva.deltaEtaJets>4.)
    {
        vbfMedium=true;
    }
    else if(vbfPreselection && mva.mDiJet>300.0 && mva.deltaEtaJets>3.)
    {
        vbfLoose=true;
    }

    bool pt0to30 = !vbfPreselection  && mva.ptDiMu <30.;
    bool pt30to50 = !vbfPreselection && mva.ptDiMu> 30. && mva.ptDiMu <50.;
    bool pt50to125 = !vbfPreselection && mva.ptDiMu> 50.  && mva.ptDiMu <125.;
    bool pt125to250 = !vbfPreselection && mva.ptDiMu> 125.  && mva.ptDiMu <250.;
    bool pt250 = !vbfPreselection && mva.ptDiMu >250.;

#ifdef BLIND
    if (!(inBlindWindow && isData))
    {
#endif
    if(!vbfPreselection)
    {
      BDTHistMuonOnly->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
      likelihoodHistMuonOnly->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      BDTHistMuonOnlyVMass->Fill(mva.mDiMu, mva.getMVA("inclusive.cfg","BDT"), weight);
      likelihoodHistMuonOnlyVMass->Fill(mva.mDiMu, mva.getMVA("inclusive.cfg","Likelihood"), weight);

      histMap2DCalibUp["BDTHistMuonOnlyVMass"]->Fill(mDiMuCalibUp, mva.getMVA("inclusive.cfg","BDT"), weight);
      histMap2DCalibDown["BDTHistMuonOnlyVMass"]->Fill(mDiMuCalibDown, mva.getMVA("inclusive.cfg","BDT"), weight);
      histMap2DResUp["BDTHistMuonOnlyVMass"]->Fill(mDiMuResUp, mva.getMVA("inclusive.cfg","BDT"), weight);
      histMap2DResDown["BDTHistMuonOnlyVMass"]->Fill(mDiMuResDown, mva.getMVA("inclusive.cfg","BDT"), weight);
      histMap2DCalibUp["likelihoodHistMuonOnlyVMass"]->Fill(mDiMuCalibUp, mva.getMVA("inclusive.cfg","Likelihood"), weight);
      histMap2DCalibDown["likelihoodHistMuonOnlyVMass"]->Fill(mDiMuCalibDown, mva.getMVA("inclusive.cfg","Likelihood"), weight);
      histMap2DResUp["likelihoodHistMuonOnlyVMass"]->Fill(mDiMuResUp, mva.getMVA("inclusive.cfg","Likelihood"), weight);
      histMap2DResDown["likelihoodHistMuonOnlyVMass"]->Fill(mDiMuResDown, mva.getMVA("inclusive.cfg","Likelihood"), weight);
    }
    else
    {
      BDTHistVBFVMass->Fill(mva.mDiMu, mva.getMVA("vbf.cfg","BDT"), weight);
      likelihoodHistVBFVMass->Fill(mva.mDiMu, mva.getMVA("vbf.cfg","Likelihood"), weight);
      BDTHistVBF->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
      likelihoodHistVBF->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);

      histMap2DCalibUp["BDTHistVBFVMass"]->Fill(mDiMuCalibUp, mva.getMVA("inclusive.cfg","BDT"), weight);
      histMap2DCalibDown["BDTHistVBFVMass"]->Fill(mDiMuCalibDown, mva.getMVA("inclusive.cfg","BDT"), weight);
      histMap2DResUp["BDTHistVBFVMass"]->Fill(mDiMuResUp, mva.getMVA("inclusive.cfg","BDT"), weight);
      histMap2DResDown["BDTHistVBFVMass"]->Fill(mDiMuResDown, mva.getMVA("inclusive.cfg","BDT"), weight);
      histMap2DCalibUp["likelihoodHistVBFVMass"]->Fill(mDiMuCalibUp, mva.getMVA("inclusive.cfg","Likelihood"), weight);
      histMap2DCalibDown["likelihoodHistVBFVMass"]->Fill(mDiMuCalibDown, mva.getMVA("inclusive.cfg","Likelihood"), weight);
      histMap2DResUp["likelihoodHistVBFVMass"]->Fill(mDiMuResUp, mva.getMVA("inclusive.cfg","Likelihood"), weight);
      histMap2DResDown["likelihoodHistVBFVMass"]->Fill(mDiMuResDown, mva.getMVA("inclusive.cfg","Likelihood"), weight);
    }
#ifdef BLIND
    }
#endif

#ifdef BLIND
    if (!(inBlindWindow && isData))
    {
#endif
    //4 GeV Window Plots
    if (mva.mDiMu < 127.0 && mva.mDiMu > 123.0)
    {
      histMap4GeVWindow["yDiMu"]->Fill(mva.yDiMu, weight);
      histMap4GeVWindow["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMap4GeVWindow["ptMu1"]->Fill(mva.ptMu1, weight);
      histMap4GeVWindow["ptMu2"]->Fill(mva.ptMu2, weight);
      histMap4GeVWindow["etaMu1"]->Fill(mva.etaMu1, weight);
      histMap4GeVWindow["etaMu2"]->Fill(mva.etaMu2, weight);
      histMap4GeVWindow["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMap4GeVWindow["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMap4GeVWindow["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMap4GeVWindow["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMap4GeVWindow["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMap4GeVWindow["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMap4GeVWindow["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMap4GeVWindow["nPU"]->Fill(nPU, weight);
      histMap4GeVWindow["nVtx"]->Fill(mva.nVtx, weight);
      histMap4GeVWindow["met"]->Fill(met.pt, weight);

      histMap4GeVWindow["mDiJet"]->Fill(mva.mDiJet, weight);
      histMap4GeVWindow["ptJet1"]->Fill(mva.ptJet1, weight);
      histMap4GeVWindow["ptJet2"]->Fill(mva.ptJet2, weight);
      histMap4GeVWindow["etaJet1"]->Fill(mva.etaJet1, weight);
      histMap4GeVWindow["etaJet2"]->Fill(mva.etaJet2, weight);
      histMap4GeVWindow["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMap4GeVWindow["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMap4GeVWindow["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMap4GeVWindow["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMap4GeVWindow["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMap4GeVWindow["nJets"]->Fill(mva.nJets, weight);
      histMap4GeVWindow["ht"]->Fill(mva.ht, weight);
      histMap4GeVWindow["mDiMu"]->Fill(mva.mDiMu, weight);
      if(!vbfPreselection)
      {
        histMap4GeVWindow["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMap4GeVWindow["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMap4GeVWindow["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMap4GeVWindow["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }
#ifdef BLIND
    }
#endif

    //DiMu Pt > 100 Plots
    if (mva.ptDiMu > 100.0)
    {
      histMapPtDiMu100["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapPtDiMu100["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapPtDiMu100["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapPtDiMu100["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapPtDiMu100["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapPtDiMu100["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapPtDiMu100["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapPtDiMu100["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapPtDiMu100["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapPtDiMu100["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapPtDiMu100["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapPtDiMu100["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapPtDiMu100["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapPtDiMu100["nPU"]->Fill(nPU, weight);
      histMapPtDiMu100["nVtx"]->Fill(mva.nVtx, weight);
      histMapPtDiMu100["met"]->Fill(met.pt, weight);

      histMapPtDiMu100["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapPtDiMu100["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapPtDiMu100["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapPtDiMu100["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapPtDiMu100["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapPtDiMu100["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapPtDiMu100["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapPtDiMu100["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapPtDiMu100["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapPtDiMu100["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapPtDiMu100["nJets"]->Fill(mva.nJets, weight);
      histMapPtDiMu100["ht"]->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapPtDiMu100["mDiMu"]->Fill(mva.mDiMu, weight);
      if(!vbfPreselection)
      {
        histMapPtDiMu100["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapPtDiMu100["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapPtDiMu100["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapPtDiMu100["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
#ifdef BLIND
      }
#endif
    }

    //VBF Preselected Plots
    if (vbfPreselection)
    {
      histMapVBFPresel["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapVBFPresel["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapVBFPresel["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapVBFPresel["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapVBFPresel["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapVBFPresel["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapVBFPresel["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapVBFPresel["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapVBFPresel["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapVBFPresel["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapVBFPresel["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapVBFPresel["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapVBFPresel["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapVBFPresel["nPU"]->Fill(nPU, weight);
      histMapVBFPresel["nVtx"]->Fill(mva.nVtx, weight);
      histMapVBFPresel["met"]->Fill(met.pt, weight);

      histMapVBFPresel["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapVBFPresel["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapVBFPresel["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapVBFPresel["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapVBFPresel["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapVBFPresel["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapVBFPresel["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapVBFPresel["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapVBFPresel["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapVBFPresel["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapVBFPresel["nJets"]->Fill(mva.nJets, weight);
      histMapVBFPresel["ht"]->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapVBFPresel["mDiMu"]->Fill(mva.mDiMu, weight);
      if(!vbfPreselection)
      {
        histMapVBFPresel["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapVBFPresel["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapVBFPresel["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapVBFPresel["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
#ifdef BLIND
      }
#endif
    }

    //Not VBF Preselected Plots
    if (!vbfPreselection)
    {
      histMapIncPresel["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapIncPresel["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapIncPresel["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapIncPresel["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapIncPresel["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapIncPresel["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapIncPresel["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapIncPresel["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapIncPresel["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapIncPresel["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapIncPresel["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapIncPresel["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapIncPresel["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapIncPresel["nPU"]->Fill(nPU, weight);
      histMapIncPresel["nVtx"]->Fill(mva.nVtx, weight);
      histMapIncPresel["met"]->Fill(met.pt, weight);

      histMapIncPresel["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapIncPresel["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapIncPresel["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapIncPresel["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapIncPresel["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapIncPresel["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapIncPresel["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapIncPresel["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapIncPresel["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapIncPresel["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapIncPresel["nJets"]->Fill(mva.nJets, weight);
      histMapIncPresel["ht"]->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapIncPresel["mDiMu"]->Fill(mva.mDiMu, weight);
      if(!vbfPreselection)
      {
        histMapIncPresel["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapIncPresel["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapIncPresel["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapIncPresel["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
#ifdef BLIND
      }
#endif
    }

    //Not In Blinding Window Plots
    if (!inBlindWindow)
    {
      histMapNotBlindWindow["mDiMu"]->Fill(mva.mDiMu, weight);
      histMapNotBlindWindow["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapNotBlindWindow["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapNotBlindWindow["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapNotBlindWindow["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapNotBlindWindow["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapNotBlindWindow["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapNotBlindWindow["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapNotBlindWindow["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapNotBlindWindow["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapNotBlindWindow["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapNotBlindWindow["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapNotBlindWindow["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapNotBlindWindow["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapNotBlindWindow["nPU"]->Fill(nPU, weight);
      histMapNotBlindWindow["nVtx"]->Fill(mva.nVtx, weight);
      histMapNotBlindWindow["met"]->Fill(met.pt, weight);

      histMapNotBlindWindow["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapNotBlindWindow["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapNotBlindWindow["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapNotBlindWindow["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapNotBlindWindow["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapNotBlindWindow["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapNotBlindWindow["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapNotBlindWindow["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapNotBlindWindow["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapNotBlindWindow["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapNotBlindWindow["nJets"]->Fill(mva.nJets, weight);
      histMapNotBlindWindow["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMapNotBlindWindow["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapNotBlindWindow["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapNotBlindWindow["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapNotBlindWindow["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }

    //Upper Control Region Plots
    if (mva.mDiMu>maxBlind)
    {
      histMapUpperControlRegion["mDiMu"]->Fill(mva.mDiMu, weight);
      histMapUpperControlRegion["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapUpperControlRegion["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapUpperControlRegion["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapUpperControlRegion["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapUpperControlRegion["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapUpperControlRegion["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapUpperControlRegion["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapUpperControlRegion["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapUpperControlRegion["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapUpperControlRegion["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapUpperControlRegion["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapUpperControlRegion["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapUpperControlRegion["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapUpperControlRegion["nPU"]->Fill(nPU, weight);
      histMapUpperControlRegion["nVtx"]->Fill(mva.nVtx, weight);
      histMapUpperControlRegion["met"]->Fill(met.pt, weight);

      histMapUpperControlRegion["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapUpperControlRegion["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapUpperControlRegion["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapUpperControlRegion["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapUpperControlRegion["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapUpperControlRegion["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapUpperControlRegion["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapUpperControlRegion["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapUpperControlRegion["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapUpperControlRegion["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapUpperControlRegion["nJets"]->Fill(mva.nJets, weight);
      histMapUpperControlRegion["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMapUpperControlRegion["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapUpperControlRegion["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapUpperControlRegion["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapUpperControlRegion["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }
    //Lower Control Region Plots
    if (mva.mDiMu<minBlind)
    {
      histMapLowerControlRegion["mDiMu"]->Fill(mva.mDiMu, weight);
      histMapLowerControlRegion["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapLowerControlRegion["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapLowerControlRegion["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapLowerControlRegion["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapLowerControlRegion["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapLowerControlRegion["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapLowerControlRegion["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapLowerControlRegion["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapLowerControlRegion["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapLowerControlRegion["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapLowerControlRegion["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapLowerControlRegion["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapLowerControlRegion["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapLowerControlRegion["nPU"]->Fill(nPU, weight);
      histMapLowerControlRegion["nVtx"]->Fill(mva.nVtx, weight);
      histMapLowerControlRegion["met"]->Fill(met.pt, weight);

      histMapLowerControlRegion["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapLowerControlRegion["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapLowerControlRegion["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapLowerControlRegion["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapLowerControlRegion["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapLowerControlRegion["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapLowerControlRegion["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapLowerControlRegion["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapLowerControlRegion["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapLowerControlRegion["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapLowerControlRegion["nJets"]->Fill(mva.nJets, weight);
      histMapLowerControlRegion["ht"]->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histMapLowerControlRegion["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
        histMapLowerControlRegion["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
      }
      else
      {
        histMapLowerControlRegion["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
        histMapLowerControlRegion["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
      }
    }

    if (vbfLoose)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapVBFLoose["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapVBFLoose["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapVBFLoose["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapVBFLoose["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapVBFLoose["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapVBFLoose["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapVBFLoose["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapVBFLoose["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapVBFLoose["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapVBFLoose["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapVBFLoose["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapVBFLoose["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapVBFLoose["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapVBFLoose["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapVBFLoose["nPU"]->Fill(nPU, weight);
      histMapVBFLoose["nVtx"]->Fill(mva.nVtx, weight);
      histMapVBFLoose["met"]->Fill(met.pt, weight);

      histMapVBFLoose["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapVBFLoose["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapVBFLoose["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapVBFLoose["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapVBFLoose["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapVBFLoose["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapVBFLoose["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapVBFLoose["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapVBFLoose["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapVBFLoose["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapVBFLoose["nJets"]->Fill(mva.nJets, weight);
      histMapVBFLoose["ht"]->Fill(mva.ht, weight);
      histMapVBFLoose["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
      histMapVBFLoose["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
    }

    if (vbfMedium)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapVBFMedium["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapVBFMedium["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapVBFMedium["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapVBFMedium["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapVBFMedium["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapVBFMedium["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapVBFMedium["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapVBFMedium["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapVBFMedium["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapVBFMedium["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapVBFMedium["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapVBFMedium["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapVBFMedium["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapVBFMedium["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapVBFMedium["nPU"]->Fill(nPU, weight);
      histMapVBFMedium["nVtx"]->Fill(mva.nVtx, weight);
      histMapVBFMedium["met"]->Fill(met.pt, weight);

      histMapVBFMedium["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapVBFMedium["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapVBFMedium["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapVBFMedium["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapVBFMedium["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapVBFMedium["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapVBFMedium["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapVBFMedium["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapVBFMedium["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapVBFMedium["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapVBFMedium["nJets"]->Fill(mva.nJets, weight);
      histMapVBFMedium["ht"]->Fill(mva.ht, weight);
      histMapVBFMedium["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
      histMapVBFMedium["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
    }

    if (vbfTight)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapVBFTight["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapVBFTight["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapVBFTight["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapVBFTight["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapVBFTight["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapVBFTight["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapVBFTight["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapVBFTight["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapVBFTight["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapVBFTight["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapVBFTight["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapVBFTight["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapVBFTight["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapVBFTight["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapVBFTight["nPU"]->Fill(nPU, weight);
      histMapVBFTight["nVtx"]->Fill(mva.nVtx, weight);
      histMapVBFTight["met"]->Fill(met.pt, weight);

      histMapVBFTight["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapVBFTight["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapVBFTight["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapVBFTight["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapVBFTight["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapVBFTight["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapVBFTight["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapVBFTight["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapVBFTight["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapVBFTight["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapVBFTight["nJets"]->Fill(mva.nJets, weight);
      histMapVBFTight["ht"]->Fill(mva.ht, weight);
      histMapVBFTight["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
      histMapVBFTight["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
    }

    if (vbfVeryTight)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapVBFVeryTight["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapVBFVeryTight["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapVBFVeryTight["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapVBFVeryTight["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapVBFVeryTight["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapVBFVeryTight["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapVBFVeryTight["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapVBFVeryTight["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapVBFVeryTight["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapVBFVeryTight["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapVBFVeryTight["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapVBFVeryTight["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapVBFVeryTight["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapVBFVeryTight["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapVBFVeryTight["nPU"]->Fill(nPU, weight);
      histMapVBFVeryTight["nVtx"]->Fill(mva.nVtx, weight);
      histMapVBFVeryTight["met"]->Fill(met.pt, weight);

      histMapVBFVeryTight["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapVBFVeryTight["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapVBFVeryTight["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapVBFVeryTight["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapVBFVeryTight["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapVBFVeryTight["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapVBFVeryTight["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapVBFVeryTight["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapVBFVeryTight["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapVBFVeryTight["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapVBFVeryTight["nJets"]->Fill(mva.nJets, weight);
      histMapVBFVeryTight["ht"]->Fill(mva.ht, weight);
      histMapVBFVeryTight["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
      histMapVBFVeryTight["likelihoodHistVBF"]->Fill(mva.getMVA("vbf.cfg","Likelihood"), weight);
    }

    if (pt0to30)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapPt0to30["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapPt0to30["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapPt0to30["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapPt0to30["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapPt0to30["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapPt0to30["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapPt0to30["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapPt0to30["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapPt0to30["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapPt0to30["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapPt0to30["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapPt0to30["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapPt0to30["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapPt0to30["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapPt0to30["nPU"]->Fill(nPU, weight);
      histMapPt0to30["nVtx"]->Fill(mva.nVtx, weight);
      histMapPt0to30["met"]->Fill(met.pt, weight);

      histMapPt0to30["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapPt0to30["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapPt0to30["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapPt0to30["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapPt0to30["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapPt0to30["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapPt0to30["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapPt0to30["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapPt0to30["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapPt0to30["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapPt0to30["nJets"]->Fill(mva.nJets, weight);
      histMapPt0to30["ht"]->Fill(mva.ht, weight);
      histMapPt0to30["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
      histMapPt0to30["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
    }

    if (pt30to50)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapPt30to50["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapPt30to50["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapPt30to50["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapPt30to50["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapPt30to50["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapPt30to50["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapPt30to50["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapPt30to50["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapPt30to50["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapPt30to50["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapPt30to50["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapPt30to50["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapPt30to50["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapPt30to50["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapPt30to50["nPU"]->Fill(nPU, weight);
      histMapPt30to50["nVtx"]->Fill(mva.nVtx, weight);
      histMapPt30to50["met"]->Fill(met.pt, weight);

      histMapPt30to50["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapPt30to50["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapPt30to50["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapPt30to50["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapPt30to50["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapPt30to50["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapPt30to50["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapPt30to50["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapPt30to50["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapPt30to50["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapPt30to50["nJets"]->Fill(mva.nJets, weight);
      histMapPt30to50["ht"]->Fill(mva.ht, weight);
      histMapPt30to50["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
      histMapPt30to50["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
    }

    if (pt50to125)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapPt50to125["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapPt50to125["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapPt50to125["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapPt50to125["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapPt50to125["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapPt50to125["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapPt50to125["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapPt50to125["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapPt50to125["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapPt50to125["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapPt50to125["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapPt50to125["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapPt50to125["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapPt50to125["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapPt50to125["nPU"]->Fill(nPU, weight);
      histMapPt50to125["nVtx"]->Fill(mva.nVtx, weight);
      histMapPt50to125["met"]->Fill(met.pt, weight);

      histMapPt50to125["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapPt50to125["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapPt50to125["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapPt50to125["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapPt50to125["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapPt50to125["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapPt50to125["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapPt50to125["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapPt50to125["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapPt50to125["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapPt50to125["nJets"]->Fill(mva.nJets, weight);
      histMapPt50to125["ht"]->Fill(mva.ht, weight);
      histMapPt50to125["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
      histMapPt50to125["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
    }

    if (pt125to250)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapPt125to250["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapPt125to250["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapPt125to250["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapPt125to250["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapPt125to250["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapPt125to250["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapPt125to250["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapPt125to250["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapPt125to250["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapPt125to250["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapPt125to250["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapPt125to250["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapPt125to250["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapPt125to250["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapPt125to250["nPU"]->Fill(nPU, weight);
      histMapPt125to250["nVtx"]->Fill(mva.nVtx, weight);
      histMapPt125to250["met"]->Fill(met.pt, weight);

      histMapPt125to250["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapPt125to250["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapPt125to250["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapPt125to250["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapPt125to250["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapPt125to250["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapPt125to250["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapPt125to250["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapPt125to250["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapPt125to250["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapPt125to250["nJets"]->Fill(mva.nJets, weight);
      histMapPt125to250["ht"]->Fill(mva.ht, weight);
      histMapPt125to250["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
      histMapPt125to250["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
    }

    if (pt250)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapPt250["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapPt250["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapPt250["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapPt250["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapPt250["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapPt250["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapPt250["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapPt250["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapPt250["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapPt250["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapPt250["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapPt250["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapPt250["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapPt250["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapPt250["nPU"]->Fill(nPU, weight);
      histMapPt250["nVtx"]->Fill(mva.nVtx, weight);
      histMapPt250["met"]->Fill(met.pt, weight);

      histMapPt250["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapPt250["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapPt250["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapPt250["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapPt250["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapPt250["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapPt250["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapPt250["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapPt250["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapPt250["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapPt250["nJets"]->Fill(mva.nJets, weight);
      histMapPt250["ht"]->Fill(mva.ht, weight);
      histMapPt250["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
      histMapPt250["likelihoodHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","Likelihood"), weight);
    }

    if (!vbfPreselection && mva.getMVAPassBDTCut("inclusive.cfg"))
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapIncBDTSig80["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapIncBDTSig80["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapIncBDTSig80["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapIncBDTSig80["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapIncBDTSig80["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapIncBDTSig80["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapIncBDTSig80["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapIncBDTSig80["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapIncBDTSig80["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapIncBDTSig80["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapIncBDTSig80["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapIncBDTSig80["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapIncBDTSig80["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapIncBDTSig80["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapIncBDTSig80["nPU"]->Fill(nPU, weight);
      histMapIncBDTSig80["nVtx"]->Fill(mva.nVtx, weight);
      histMapIncBDTSig80["met"]->Fill(met.pt, weight);

      histMapIncBDTSig80["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapIncBDTSig80["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapIncBDTSig80["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapIncBDTSig80["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapIncBDTSig80["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapIncBDTSig80["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapIncBDTSig80["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapIncBDTSig80["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapIncBDTSig80["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapIncBDTSig80["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapIncBDTSig80["nJets"]->Fill(mva.nJets, weight);
      histMapIncBDTSig80["ht"]->Fill(mva.ht, weight);
      histMapIncBDTSig80["BDTHistMuonOnly"]->Fill(mva.getMVA("inclusive.cfg","BDT"), weight);
    }

    if (vbfPreselection && mva.getMVAPassBDTCut("vbf.cfg"))
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histMapVBFBDTSig80["mDiMu"]->Fill(mva.mDiMu, weight);
#ifdef BLIND
      }
#endif
      histMapVBFBDTSig80["yDiMu"]->Fill(mva.yDiMu, weight);
      histMapVBFBDTSig80["ptDiMu"]->Fill(mva.ptDiMu, weight);
      histMapVBFBDTSig80["ptMu1"]->Fill(mva.ptMu1, weight);
      histMapVBFBDTSig80["ptMu2"]->Fill(mva.ptMu2, weight);
      histMapVBFBDTSig80["etaMu1"]->Fill(mva.etaMu1, weight);
      histMapVBFBDTSig80["etaMu2"]->Fill(mva.etaMu2, weight);
      histMapVBFBDTSig80["cosThetaStar"]->Fill(mva.cosThetaStar, weight);
      histMapVBFBDTSig80["cosThetaStarCS"]->Fill(mva.cosThetaStarCS, weight);
      histMapVBFBDTSig80["deltaPhiMuons"]->Fill(mva.deltaPhiMuons, weight);
      histMapVBFBDTSig80["deltaEtaMuons"]->Fill(mva.deltaEtaMuons, weight);
      histMapVBFBDTSig80["deltaRMuons"]->Fill(mva.deltaRMuons, weight);
      histMapVBFBDTSig80["relIsoMu1"]->Fill(mva.relIsoMu1, weight);
      histMapVBFBDTSig80["relIsoMu2"]->Fill(mva.relIsoMu2, weight);
      histMapVBFBDTSig80["nPU"]->Fill(nPU, weight);
      histMapVBFBDTSig80["nVtx"]->Fill(mva.nVtx, weight);
      histMapVBFBDTSig80["met"]->Fill(met.pt, weight);

      histMapVBFBDTSig80["mDiJet"]->Fill(mva.mDiJet, weight);
      histMapVBFBDTSig80["ptJet1"]->Fill(mva.ptJet1, weight);
      histMapVBFBDTSig80["ptJet2"]->Fill(mva.ptJet2, weight);
      histMapVBFBDTSig80["etaJet1"]->Fill(mva.etaJet1, weight);
      histMapVBFBDTSig80["etaJet2"]->Fill(mva.etaJet2, weight);
      histMapVBFBDTSig80["deltaEtaJets"]->Fill(mva.deltaEtaJets, weight);
      histMapVBFBDTSig80["deltaPhiJets"]->Fill(mva.deltaPhiJets, weight);
      histMapVBFBDTSig80["deltaRJets"]->Fill(mva.deltaRJets, weight);
      histMapVBFBDTSig80["nJetsInRapidityGap"]->Fill(mva.nJetsInRapidityGap, weight);
      histMapVBFBDTSig80["htInRapidityGap"]->Fill(mva.htInRapidityGap, weight);
      histMapVBFBDTSig80["nJets"]->Fill(mva.nJets, weight);
      histMapVBFBDTSig80["ht"]->Fill(mva.ht, weight);
      histMapVBFBDTSig80["BDTHistVBF"]->Fill(mva.getMVA("vbf.cfg","BDT"), weight);
    }

  }// end event loop

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  outFile->cd();

  for(histMapIter = histMap.begin(); histMapIter != histMap.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  for(histMap2DIter = histMap2D.begin(); histMap2DIter != histMap2D.end(); histMap2DIter++)
  {
    histMap2DIter->second->Write();
  }

  TDirectory* dirCalibUp = outFile->mkdir("CalibUp");
  dirCalibUp->cd();
  for(histMap2DIter = histMap2DCalibUp.begin(); histMap2DIter != histMap2DCalibUp.end(); histMap2DIter++)
  {
    histMap2DIter->second->Write();
  }

  TDirectory* dirCalibDown = outFile->mkdir("CalibDown");
  dirCalibDown->cd();
  for(histMap2DIter = histMap2DCalibDown.begin(); histMap2DIter != histMap2DCalibDown.end(); histMap2DIter++)
  {
    histMap2DIter->second->Write();
  }

  TDirectory* dirResUp = outFile->mkdir("ResUp");
  dirResUp->cd();
  for(histMap2DIter = histMap2DResUp.begin(); histMap2DIter != histMap2DResUp.end(); histMap2DIter++)
  {
    histMap2DIter->second->Write();
  }

  TDirectory* dirResDown = outFile->mkdir("ResDown");
  dirResDown->cd();
  for(histMap2DIter = histMap2DResDown.begin(); histMap2DIter != histMap2DResDown.end(); histMap2DIter++)
  {
    histMap2DIter->second->Write();
  }

  TDirectory* dir4GeVWindow = outFile->mkdir("4GeVWindow");
  dir4GeVWindow->cd();
  for(histMapIter = histMap4GeVWindow.begin(); histMapIter != histMap4GeVWindow.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirPtDiMu100 = outFile->mkdir("PtDiMu100");
  dirPtDiMu100->cd();
  for(histMapIter = histMapPtDiMu100.begin(); histMapIter != histMapPtDiMu100.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirVBFPresel = outFile->mkdir("VBFPresel");
  dirVBFPresel->cd();
  for(histMapIter = histMapVBFPresel.begin(); histMapIter != histMapVBFPresel.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirIncPresel = outFile->mkdir("IncPresel");
  dirIncPresel->cd();
  for(histMapIter = histMapIncPresel.begin(); histMapIter != histMapIncPresel.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirNotBlindWindow = outFile->mkdir("NotBlindWindow");
  dirNotBlindWindow->cd();
  for(histMapIter = histMapNotBlindWindow.begin(); histMapIter != histMapNotBlindWindow.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirUpperControlRegion = outFile->mkdir("UpperControlRegion");
  dirUpperControlRegion->cd();
  for(histMapIter = histMapUpperControlRegion.begin(); histMapIter != histMapUpperControlRegion.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirLowerControlRegion = outFile->mkdir("LowerControlRegion");
  dirLowerControlRegion->cd();
  for(histMapIter = histMapLowerControlRegion.begin(); histMapIter != histMapLowerControlRegion.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirVBFLoose = outFile->mkdir("VBFLoose");
  dirVBFLoose->cd();
  for(histMapIter = histMapVBFLoose.begin(); histMapIter != histMapVBFLoose.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirVBFMedium = outFile->mkdir("VBFMedium");
  dirVBFMedium->cd();
  for(histMapIter = histMapVBFMedium.begin(); histMapIter != histMapVBFMedium.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirVBFTight = outFile->mkdir("VBFTight");
  dirVBFTight->cd();
  for(histMapIter = histMapVBFTight.begin(); histMapIter != histMapVBFTight.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirVBFVeryTight = outFile->mkdir("VBFVeryTight");
  dirVBFVeryTight->cd();
  for(histMapIter = histMapVBFVeryTight.begin(); histMapIter != histMapVBFVeryTight.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirPt0to30 = outFile->mkdir("Pt0to30");
  dirPt0to30->cd();
  for(histMapIter = histMapPt0to30.begin(); histMapIter != histMapPt0to30.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirPt30to50 = outFile->mkdir("Pt30to50");
  dirPt30to50->cd();
  for(histMapIter = histMapPt30to50.begin(); histMapIter != histMapPt30to50.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirPt50to125 = outFile->mkdir("Pt50to125");
  dirPt50to125->cd();
  for(histMapIter = histMapPt50to125.begin(); histMapIter != histMapPt50to125.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirPt125to250 = outFile->mkdir("Pt125to250");
  dirPt125to250->cd();
  for(histMapIter = histMapPt125to250.begin(); histMapIter != histMapPt125to250.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirPt250 = outFile->mkdir("Pt250");
  dirPt250->cd();
  for(histMapIter = histMapPt250.begin(); histMapIter != histMapPt250.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirIncBDTSig80 = outFile->mkdir("IncBDTSig80");
  dirIncBDTSig80->cd();
  for(histMapIter = histMapIncBDTSig80.begin(); histMapIter != histMapIncBDTSig80.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  TDirectory* dirVBFBDTSig80 = outFile->mkdir("VBFBDTSig80");
  dirVBFBDTSig80->cd();
  for(histMapIter = histMapVBFBDTSig80.begin(); histMapIter != histMapVBFBDTSig80.end(); histMapIter++)
  {
    histMapIter->second->Write();
  }

  ofstream testOutFile;
  testOutFile.open("testEventNums.txt");
  testOutFile << testString;
  testOutFile.close();
  cout <<"#######################\n"<< testString <<"#######################\n" << endl;
  cout << "testCounter: "<< testCounter << endl;
  cout << "analyzer done." << endl << endl;
  return 0;
}

void
fillMuonHist(TH1F* hist, _MuonInfo& mu1, _MuonInfo& mu2)
{

  hist->Fill(1.0);
  if (!mu1.isGlobal || !mu2.isGlobal)  return;
  hist->Fill(2.0);
  if (!mu1.isPFMuon || !mu2.isPFMuon) return;
  hist->Fill(3.0);

  // acceptance cuts
  if (mu1.pt < 25 || mu2.pt < 25)         return; // pt cut
  hist->Fill(4.0);
  if (fabs(mu1.eta) > 2.1 || fabs(mu2.eta) > 2.1) return; // eta cut
  hist->Fill(5.0);

  // kinematic cuts
  if (mu1.numTrackerLayers < 6 || mu2.numTrackerLayers < 6) return; // # hits in tracker
  hist->Fill(6.0);

  if(getRelIso(mu1) > 0.12 || getRelIso(mu2) > 0.12)
  {
      //cout << "Iso 1: "<< getRelIso(mu1) << "    Iso 2: " << getRelIso(mu2) << endl;
      return;
  }
  hist->Fill(7.0);

  if (fabs(mu1.d0_PV) > 0.2 || fabs(mu2.d0_PV) > 0.2) return;
  hist->Fill(8.0);
  if (fabs(mu1.dz_PV) > 0.5 || fabs(mu2.dz_PV) > 0.5) return;
  hist->Fill(9.0);

  if ( mu1.numValidMuonHits  < 1  || mu2.numValidMuonHits  < 1) return;
  hist->Fill(10.0);
  if ( mu1.numValidPixelHits < 1  || mu2.numValidPixelHits < 1) return;
  hist->Fill(11.0);
  if ( mu1.numOfMatchedStations < 2  || mu2.numOfMatchedStations < 2)
  {
      //cout << "Sta 1: "<<mu1.numOfMatchedStations << "    Sta 2: " << mu2.numOfMatchedStations << endl;
      return;
  }
  hist->Fill(12.0);
  if ( mu1.normChiSquare > 10 || mu2.normChiSquare > 10)     return;
  hist->Fill(13.0);

}

void
printStationMiss(_MuonInfo& mu1, _MuonInfo& mu2, _EventInfo& eventInfo, std::string & testString, unsigned & testCounter)
{

  if (!mu1.isGlobal || !mu2.isGlobal)  return;
  if (!mu1.isPFMuon || !mu2.isPFMuon) return;

  // acceptance cuts
  if (mu1.pt < 25 || mu2.pt < 25)         return; // pt cut
  if (fabs(mu1.eta) > 2.1 || fabs(mu2.eta) > 2.1) return; // eta cut

  // kinematic cuts
  if (mu1.numTrackerLayers < 6 || mu2.numTrackerLayers < 6) return; // # hits in tracker

  if(getRelIso(mu1) > 0.12 || getRelIso(mu2) > 0.12)
  {
      //cout << "Iso 1: "<< getRelIso(mu1) << "    Iso 2: " << getRelIso(mu2) << endl;
      return;
  }

  if (fabs(mu1.d0_PV) > 0.2 || fabs(mu2.d0_PV) > 0.2) return;
  if (fabs(mu1.dz_PV) > 0.5 || fabs(mu2.dz_PV) > 0.5) return;

  if ( mu1.numValidMuonHits  < 1  || mu2.numValidMuonHits  < 1) return;
  if ( mu1.numValidPixelHits < 1  || mu2.numValidPixelHits < 1) return;
  if ( mu1.numOfMatchedStations < 2  || mu2.numOfMatchedStations < 2)
  {
        testCounter++;
        //std::cout <<eventInfo.run <<":"<<eventInfo.event <<"\n"<< std::endl;
        testString.appendAny(eventInfo.run);
        testString.append(":");
        testString.appendAny(eventInfo.event);
        //testString.append("  #  ");
        //testString.appendAny(mu1.numOfMatchedStations);
        //testString.append("  ");
        //testString.appendAny(mu2.numOfMatchedStations);
        testString.append("\n");
      
      return;
  }
  if ( mu1.normChiSquare > 10 || mu2.normChiSquare > 10)     return;

}
