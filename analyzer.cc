#include <algorithm>
#include <limits.h>
#include <ctime>
#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include <boost/lexical_cast.hpp>
//Defines method of std::string that appends any type :-)
#define appendAny(a) append(boost::lexical_cast<std::string>(a))

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

#define JETPUID
#define PUREWEIGHT
//#define SMEARING

using namespace std;
using namespace boost;

void fillMuonHist(TH1F* hist, _MuonInfo& mu1, _MuonInfo& mu2);
void printStationMiss(_MuonInfo& mu1, _MuonInfo& mu2, _EventInfo& eventInfo, std::string & testString, unsigned & testCounter);

struct HistStruct
{
  HistStruct();
  ~HistStruct();
  void Write();

  std::vector<TH1F*> histVec;
  std::vector<TH2F*> histVec2D;

  TH1F* mDiMu;
  TH1F* mDiMuBB;
  TH1F* mDiMuBO;
  TH1F* mDiMuBE;
  TH1F* mDiMuOO;
  TH1F* mDiMuOE;
  TH1F* mDiMuEE;
  TH1F* mDiMuNotBB;

  TH1F* mDiJet;
  TH1F* ptDiMu;
  TH1F* ptDiJet;
  TH1F* yDiMu;

  TH2F* yVptDiMu;
  TH2F* ptVmDiMu;
  TH2F* yVmDiMu;
  TH2F* phiVmDiMu;

  TH1F* ptMu1;
  TH1F* ptMu2;
  TH1F* ptJet1;
  TH1F* ptJet2;

  TH1F* etaMu1;
  TH1F* etaMu2;
  TH1F* etaJet1;
  TH1F* etaJet2;

  TH1F* deltaEtaJets;
  TH1F* deltaPhiJets;
  TH1F* deltaRJets;

  TH1F* deltaEtaMuons;
  TH1F* deltaPhiMuons;
  TH1F* deltaRMuons;

  TH1F* countsHist;
  TH1F* countsHist2;

  TH1F* cosThetaStar;
  TH1F* cosThetaStarCS;

  TH1F* puJetIDSimpleDiscJet1;
  TH1F* puJetIDSimpleDiscJet2;
  TH1F* puJetIDSimpleDiscJet3;

  TH1F* puJetIDSimpleJet1;
  TH1F* puJetIDSimpleJet2;
  TH1F* puJetIDSimpleJet3;

  TH1F* BDTHistMuonOnly;
  TH1F* likelihoodHistMuonOnly;

  TH1F* BDTHistVBF;
  TH1F* likelihoodHistVBF;

  TH2F* BDTHistMuonOnlyVMass;
  TH2F* likelihoodHistMuonOnlyVMass;
  TH2F* BDTHistVBFVMass;
  TH2F* likelihoodHistVBFVMass;

  TH1F* relIsoMu1;
  TH1F* relIsoMu2;

  TH1F* nJets;
  TH1F* ht;
  TH1F* nJetsInRapidityGap;
  TH1F* htInRapidityGap;

  TH1F* nPU;
  TH1F* nVtx;
  TH1F* met;
  TH1F* weight;
};

int main(int argc, char *argv[])
{
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  gErrorIgnoreLevel = kError;
  time_t timeStart = time(NULL);

  const char* optionIntro = "H->MuMu Analyzer\n\nUsage: ./analyzer [--help] [--train] [--maxEvents N] <outputFileName.root> <inputFileName.root> [<inputFileName2.root>...]\n\nAllowed Options";
  program_options::options_description optionDesc(optionIntro);
  optionDesc.add_options()
      ("help,h", "Produce Help Message")
      ("trainingTree,t",program_options::value<string>(), "Create Training Tree File with filename")
      ("filenames",program_options::value<vector<string> >(), "Input & Output File Names, put output name first followed by all input file names")
      ("maxEvents,m",program_options::value<int>(), "Maximum Number of Events to Process")
      ("runPeriod,r",program_options::value<string>(), "Running Perios e.g. 7TeV, 8TeV")
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

  string runPeriod = "";
  if (optionMap.count("runPeriod"))
  {
    runPeriod = optionMap["runPeriod"].as<string>();
  }
  //else
  //{
  //  cout << "Run Period not specified, exiting"<< endl;
  //  return 1;
  //}
  cout << "Run Period: " << runPeriod << endl;

  /////////////////////////////
  //////////// Setup //////////
  /////////////////////////////

  float minMmm = 70.0;
  float maxMmm = 160.0;
  //float minMmm = 110.0;
  //float maxMmm = 150.0;

  float minBlind = 120;
  float maxBlind = 130;

  float calib = -0.1;
  float calibSysSmear = 0.2;
  float resSmear = 1.169; // should all be around 1; a ratio of resolutions
  float resSysSmear = 0.2; // Error on that ratio

  string cfgNameInc = "inclusive_";
  string cfgNameVBF = "vbf_";
  cfgNameInc.append(runPeriod);
  cfgNameVBF.append(runPeriod);
  cfgNameInc.append(".cfg");
  cfgNameVBF.append(".cfg");
  cout << "cfgNames: " << cfgNameInc <<" \t"<< cfgNameVBF << endl;

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

  ///////////////////////////////
  // Which Muon Selection to Use

  bool (*muonIdFuncPtr)(_MuonInfo&);
  if (runPeriod == "7TeV")
  {
    cout << "Using 2011 Tight Muon Selection\n";
    muonIdFuncPtr = &isKinTight_2011;
  }
  else
  {
    cout << "Using 2012 Tight Muon Selection\n";
    muonIdFuncPtr = &isKinTight_2012;
  }

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
  bool trainingTreeRun=false;
  if (optionMap.count("trainingTree")) 
  {
      cout << "Training enabled" << "\n";
      trainingTreeFileName = optionMap["trainingTree"].as<string>();
      cout << "Training Tree File Name: " << trainingTreeFileName << "\n";
      trainingTreeRun=true;
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
  HistStruct hists;
  HistStruct hists4GeVWindow;
  HistStruct histsPtDiMu100;
  HistStruct histsVBFPresel;
  HistStruct histsIncPresel;
  HistStruct histsNotBlindWindow;
  HistStruct histsUpperControlRegion;
  HistStruct histsLowerControlRegion;
  HistStruct histsVBFLoose;
  HistStruct histsVBFMedium;
  HistStruct histsVBFTight;
  HistStruct histsVBFVeryTight;
  HistStruct histsPt0to30;
  HistStruct histsPt30to50;
  HistStruct histsPt50to125;
  HistStruct histsPt125to250;
  HistStruct histsPt250;
  HistStruct histsIncBDTSig80;
  HistStruct histsVBFBDTSig80;
  HistStruct histsCalibUp;
  HistStruct histsCalibDown;
  HistStruct histsResUp;
  HistStruct histsResDown;

  //////////////////////////
  //for MVA

  std::vector<std::string> mvaConfigNames;
  mvaConfigNames.push_back(cfgNameInc);
  mvaConfigNames.push_back(cfgNameVBF);
  MVA mva(mvaConfigNames,trainingTreeFileName);

  TRandom3 random(1457);

  //////////////////////////
  //for PU reweighting

#ifdef PUREWEIGHT
  reweight::LumiReWeighting lumiWeights("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root","pileupDists/PileUpHist2012ABC.root","pileup","pileup");
  if (runPeriod == "7TeV")
  {
    cout << "Using 2011AB PU reweighting\n";
    lumiWeights = reweight::LumiReWeighting("pileupDists/PileUpHistMCFall11.root","pileupDists/PileUpHist2011AB.root","pileup","pileup");
  }
  else
  {
    cout << "Using 2012ABC PU reweighting\n";
  }
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
  
  float timeReading = 0.;
  float timeReadingAll = 0.;
  float timeProcessing = 0.;
  float timeFilling = 0.;
  time_t timeStartEventLoop = time(NULL);
  for(unsigned i=0; i<nEvents;i++)
  {
    if(i >= maxEvents)
      continue;
    time_t timeStartReading = time(NULL);
    tree->GetEvent(i);
    time_t timeStopReading = time(NULL);
    timeReadingAll += difftime(timeStopReading,timeStartReading);
    if (i % reportEach == 0) cout << "Event: " << i << endl;

    fillMuonHist(hists.countsHist2, reco1, reco2);
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

    hists.countsHist->Fill(0.0, weight);

    if (!((*muonIdFuncPtr)(reco1)) || !((*muonIdFuncPtr)(reco2)))
          continue;

    hists.countsHist->Fill(1.0, weight);

    if (!isHltMatched(reco1,reco2,allowedHLTPaths) && isData)
        continue;

    hists.countsHist->Fill(2.0, weight);

    if (reco1.charge*reco2.charge != -1)
        continue;

    hists.countsHist->Fill(3.0, weight);

    if (mva.mDiMu < minMmm || mva.mDiMu > maxMmm)
        continue;

    hists.countsHist->Fill(4.0, weight);

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

    bool isBB = false;
    bool isBO = false;
    bool isBE = false;
    bool isOO = false;
    bool isOE = false;
    bool isEE = false;
    if(fabs(muon1.eta)<0.8 && fabs(muon2.eta)<0.8)
    {
        isBB=true;
    }
    else if(
        (fabs(muon1.eta)<0.8 && fabs(muon2.eta)<1.6)
            || (fabs(muon1.eta)<1.6 && fabs(muon2.eta)<0.8)
        )
    {
        isBO=true;
    }
    else if(
        fabs(muon1.eta)<0.8 || fabs(muon2.eta)<0.8
        )
    {
        isBE=true;
    }
    else if(
        fabs(muon1.eta)<1.6 && fabs(muon2.eta)<1.6
        )
    {
        isOO=true;
    }
    else if(
        fabs(muon1.eta)<1.6 || fabs(muon2.eta)<1.6
        )
    {
        isOE=true;
    }
    else
    {
        isEE=true;
    }
    bool isNotBB = !isBB;
    

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
    hists.mDiMu->Fill(mva.mDiMu, weight);
    if (isBB) hists.mDiMuBB->Fill(mva.mDiMu, weight);
    if (isBE) hists.mDiMuBE->Fill(mva.mDiMu, weight);
    if (isBO) hists.mDiMuBO->Fill(mva.mDiMu, weight);
    if (isOO) hists.mDiMuOO->Fill(mva.mDiMu, weight);
    if (isOE) hists.mDiMuOE->Fill(mva.mDiMu, weight);
    if (isEE) hists.mDiMuEE->Fill(mva.mDiMu, weight);
    if (isNotBB) hists.mDiMuNotBB->Fill(mva.mDiMu, weight);

    hists.yVmDiMu->Fill(mva.mDiMu,fabs(mva.yDiMu), weight);
    hists.ptVmDiMu->Fill(mva.mDiMu,mva.ptDiMu, weight);
    hists.phiVmDiMu->Fill(mva.mDiMu,recoCandPhi, weight);
#ifdef BLIND
    }
#endif

    hists.yDiMu->Fill(mva.yDiMu, weight);
    hists.ptDiMu->Fill(mva.ptDiMu, weight);
    hists.yVptDiMu->Fill(mva.ptDiMu,fabs(mva.yDiMu), weight);
    hists.ptMu1->Fill(mva.ptMu1, weight);
    hists.ptMu2->Fill(mva.ptMu2, weight);
    hists.etaMu1->Fill(mva.etaMu1, weight);
    hists.etaMu2->Fill(mva.etaMu2, weight);
    hists.cosThetaStar->Fill(mva.cosThetaStar, weight);
    hists.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);

    hists.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
    hists.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
    hists.deltaRMuons->Fill(mva.deltaRMuons, weight);

    hists.relIsoMu1->Fill(mva.relIsoMu1, weight);
    hists.relIsoMu2->Fill(mva.relIsoMu2, weight);

    hists.nPU->Fill(nPU, weight);
    hists.nVtx->Fill(mva.nVtx, weight);
    hists.met->Fill(met.pt, weight);
    hists.weight->Fill(weight);

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
            hists.puJetIDSimpleDiscJet1->Fill(puJetSimpleDisc[iJet], weight);
          else if (iJet==1)
            hists.puJetIDSimpleDiscJet2->Fill(puJetSimpleDisc[iJet], weight);
          else if (iJet==2)
            hists.puJetIDSimpleDiscJet3->Fill(puJetSimpleDisc[iJet], weight);

          if (iJet==0)
            hists.puJetIDSimpleJet1->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose), weight);
          else if (iJet==1)
            hists.puJetIDSimpleJet2->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose), weight);
          else if (iJet==2)
            hists.puJetIDSimpleJet3->Fill(passPUJetID(puJetSimpleId[iJet],puJetLoose), weight);
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

      hists.mDiJet->Fill(mva.mDiJet, weight);
      hists.ptDiJet->Fill(mva.ptDiJet, weight);
      hists.ptJet1->Fill(mva.ptJet1, weight);
      hists.ptJet2->Fill(mva.ptJet2, weight);
      hists.etaJet1->Fill(mva.etaJet1, weight);
      hists.etaJet2->Fill(mva.etaJet2, weight);
      hists.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      hists.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      hists.deltaRJets->Fill(mva.deltaRJets, weight);
      hists.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      hists.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      hists.nJets->Fill(mva.nJets, weight);
      hists.ht->Fill(mva.ht, weight);
    }
  
//HIG-12-007 PAS H->tautau
//The VBF category requires at least two jets with pT > 30 GeV/c, |η1 − η2 | > 4.0 and
//η1 · η2 < 0 (with η1 the pseudorapidty of the leading jet and η2 the pseudorapidity
//of the subleading jet), and a di-jet invariant mass m12 > 400 GeV/c2 , with no other
//jet with pT > 30 GeV/c in the rapidity region between the two jets.

    mva.writeEvent();

    if (trainingTreeRun) //Skip Filling of histos when training Tree
        continue;

    bool vbfPreselection = mva.mDiJet>300.0 && mva.deltaEtaJets>3.0 && mva.productEtaJets<0.0 && mva.nJetsInRapidityGap == 0;
    //if(vbfPreselection)
    //  std::cout << "VBF Preselected!!";

    if(vbfPreselection)
        hists.countsHist->Fill(6);
    else
        hists.countsHist->Fill(5);
    

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

    float bdtValInc =  mva.getMVA(cfgNameInc,"BDT");
    float bdtValVBF =  mva.getMVA(cfgNameVBF,"BDT");
    float likeValInc =  mva.getMVA(cfgNameInc,"Likelihood");
    float likeValVBF =  mva.getMVA(cfgNameVBF,"Likelihood");

    time_t timeStartFilling = time(NULL);

#ifdef BLIND
    if (!(inBlindWindow && isData))
    {
#endif
    if(!vbfPreselection)
    {
      hists.BDTHistMuonOnly->Fill(bdtValInc, weight);
      hists.likelihoodHistMuonOnly->Fill(likeValInc, weight);
      hists.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      hists.likelihoodHistMuonOnlyVMass->Fill(mva.mDiMu, likeValInc, weight);

      histsCalibUp.BDTHistMuonOnlyVMass->Fill(mDiMuCalibUp, bdtValInc, weight);
      histsCalibDown.BDTHistMuonOnlyVMass->Fill(mDiMuCalibDown, bdtValInc, weight);
      histsResUp.BDTHistMuonOnlyVMass->Fill(mDiMuResUp, bdtValInc, weight);
      histsResDown.BDTHistMuonOnlyVMass->Fill(mDiMuResDown, bdtValInc, weight);
      histsCalibUp.likelihoodHistMuonOnlyVMass->Fill(mDiMuCalibUp, likeValInc, weight);
      histsCalibDown.likelihoodHistMuonOnlyVMass->Fill(mDiMuCalibDown, likeValInc, weight);
      histsResUp.likelihoodHistMuonOnlyVMass->Fill(mDiMuResUp, likeValInc, weight);
      histsResDown.likelihoodHistMuonOnlyVMass->Fill(mDiMuResDown, likeValInc, weight);
    }
    else
    {
      hists.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      hists.likelihoodHistVBFVMass->Fill(mva.mDiMu, likeValVBF, weight);
      hists.BDTHistVBF->Fill(bdtValVBF, weight);
      hists.likelihoodHistVBF->Fill(likeValVBF, weight);

      histsCalibUp.BDTHistVBFVMass->Fill(mDiMuCalibUp, bdtValInc, weight);
      histsCalibDown.BDTHistVBFVMass->Fill(mDiMuCalibDown, bdtValInc, weight);
      histsResUp.BDTHistVBFVMass->Fill(mDiMuResUp, bdtValInc, weight);
      histsResDown.BDTHistVBFVMass->Fill(mDiMuResDown, bdtValInc, weight);
      histsCalibUp.likelihoodHistVBFVMass->Fill(mDiMuCalibUp, likeValInc, weight);
      histsCalibDown.likelihoodHistVBFVMass->Fill(mDiMuCalibDown, likeValInc, weight);
      histsResUp.likelihoodHistVBFVMass->Fill(mDiMuResUp, likeValInc, weight);
      histsResDown.likelihoodHistVBFVMass->Fill(mDiMuResDown, likeValInc, weight);
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
      hists4GeVWindow.yDiMu->Fill(mva.yDiMu, weight);
      hists4GeVWindow.ptDiMu->Fill(mva.ptDiMu, weight);
      hists4GeVWindow.ptMu1->Fill(mva.ptMu1, weight);
      hists4GeVWindow.ptMu2->Fill(mva.ptMu2, weight);
      hists4GeVWindow.etaMu1->Fill(mva.etaMu1, weight);
      hists4GeVWindow.etaMu2->Fill(mva.etaMu2, weight);
      hists4GeVWindow.cosThetaStar->Fill(mva.cosThetaStar, weight);
      hists4GeVWindow.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      hists4GeVWindow.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      hists4GeVWindow.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      hists4GeVWindow.deltaRMuons->Fill(mva.deltaRMuons, weight);
      hists4GeVWindow.relIsoMu1->Fill(mva.relIsoMu1, weight);
      hists4GeVWindow.relIsoMu2->Fill(mva.relIsoMu2, weight);
      hists4GeVWindow.nPU->Fill(nPU, weight);
      hists4GeVWindow.nVtx->Fill(mva.nVtx, weight);
      hists4GeVWindow.met->Fill(met.pt, weight);

      hists4GeVWindow.mDiJet->Fill(mva.mDiJet, weight);
      hists4GeVWindow.ptDiJet->Fill(mva.ptDiJet, weight);
      hists4GeVWindow.ptJet1->Fill(mva.ptJet1, weight);
      hists4GeVWindow.ptJet2->Fill(mva.ptJet2, weight);
      hists4GeVWindow.etaJet1->Fill(mva.etaJet1, weight);
      hists4GeVWindow.etaJet2->Fill(mva.etaJet2, weight);
      hists4GeVWindow.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      hists4GeVWindow.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      hists4GeVWindow.deltaRJets->Fill(mva.deltaRJets, weight);
      hists4GeVWindow.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      hists4GeVWindow.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      hists4GeVWindow.nJets->Fill(mva.nJets, weight);
      hists4GeVWindow.ht->Fill(mva.ht, weight);
      hists4GeVWindow.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) hists4GeVWindow.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) hists4GeVWindow.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) hists4GeVWindow.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) hists4GeVWindow.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) hists4GeVWindow.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) hists4GeVWindow.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) hists4GeVWindow.mDiMuNotBB->Fill(mva.mDiMu, weight);

      if(!vbfPreselection)
      {
        hists4GeVWindow.BDTHistMuonOnly->Fill(bdtValInc, weight);
        hists4GeVWindow.likelihoodHistMuonOnly->Fill(likeValInc, weight);
      }
      else
      {
        hists4GeVWindow.BDTHistVBF->Fill(bdtValVBF, weight);
        hists4GeVWindow.likelihoodHistVBF->Fill(likeValVBF, weight);
      }
    }
#ifdef BLIND
    }
#endif

    //DiMu Pt > 100 Plots
    if (mva.ptDiMu > 100.0)
    {
      histsPtDiMu100.yDiMu->Fill(mva.yDiMu, weight);
      histsPtDiMu100.ptDiMu->Fill(mva.ptDiMu, weight);
      histsPtDiMu100.ptMu1->Fill(mva.ptMu1, weight);
      histsPtDiMu100.ptMu2->Fill(mva.ptMu2, weight);
      histsPtDiMu100.etaMu1->Fill(mva.etaMu1, weight);
      histsPtDiMu100.etaMu2->Fill(mva.etaMu2, weight);
      histsPtDiMu100.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsPtDiMu100.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsPtDiMu100.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsPtDiMu100.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsPtDiMu100.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsPtDiMu100.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsPtDiMu100.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsPtDiMu100.nPU->Fill(nPU, weight);
      histsPtDiMu100.nVtx->Fill(mva.nVtx, weight);
      histsPtDiMu100.met->Fill(met.pt, weight);

      histsPtDiMu100.mDiJet->Fill(mva.mDiJet, weight);
      histsPtDiMu100.ptDiJet->Fill(mva.ptDiJet, weight);
      histsPtDiMu100.ptJet1->Fill(mva.ptJet1, weight);
      histsPtDiMu100.ptJet2->Fill(mva.ptJet2, weight);
      histsPtDiMu100.etaJet1->Fill(mva.etaJet1, weight);
      histsPtDiMu100.etaJet2->Fill(mva.etaJet2, weight);
      histsPtDiMu100.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsPtDiMu100.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsPtDiMu100.deltaRJets->Fill(mva.deltaRJets, weight);
      histsPtDiMu100.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsPtDiMu100.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsPtDiMu100.nJets->Fill(mva.nJets, weight);
      histsPtDiMu100.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsPtDiMu100.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsPtDiMu100.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsPtDiMu100.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsPtDiMu100.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsPtDiMu100.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsPtDiMu100.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsPtDiMu100.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsPtDiMu100.mDiMuNotBB->Fill(mva.mDiMu, weight);

      if(!vbfPreselection)
      {
        histsPtDiMu100.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsPtDiMu100.likelihoodHistMuonOnly->Fill(likeValInc, weight);
      }
      else
      {
        histsPtDiMu100.BDTHistVBF->Fill(bdtValVBF, weight);
        histsPtDiMu100.likelihoodHistVBF->Fill(likeValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    //VBF Preselected Plots
    if (vbfPreselection)
    {
      histsVBFPresel.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFPresel.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFPresel.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFPresel.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFPresel.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFPresel.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFPresel.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFPresel.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFPresel.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFPresel.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFPresel.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFPresel.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFPresel.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFPresel.nPU->Fill(nPU, weight);
      histsVBFPresel.nVtx->Fill(mva.nVtx, weight);
      histsVBFPresel.met->Fill(met.pt, weight);

      histsVBFPresel.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFPresel.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFPresel.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFPresel.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFPresel.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFPresel.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFPresel.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFPresel.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFPresel.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFPresel.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFPresel.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFPresel.nJets->Fill(mva.nJets, weight);
      histsVBFPresel.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFPresel.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsVBFPresel.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsVBFPresel.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsVBFPresel.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsVBFPresel.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsVBFPresel.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsVBFPresel.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsVBFPresel.mDiMuNotBB->Fill(mva.mDiMu, weight);

      if(!vbfPreselection)
      {
        histsVBFPresel.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsVBFPresel.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsVBFPresel.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsVBFPresel.BDTHistVBF->Fill(bdtValVBF, weight);
        histsVBFPresel.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsVBFPresel.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    //Not VBF Preselected Plots
    if (!vbfPreselection)
    {
      histsIncPresel.yDiMu->Fill(mva.yDiMu, weight);
      histsIncPresel.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncPresel.ptMu1->Fill(mva.ptMu1, weight);
      histsIncPresel.ptMu2->Fill(mva.ptMu2, weight);
      histsIncPresel.etaMu1->Fill(mva.etaMu1, weight);
      histsIncPresel.etaMu2->Fill(mva.etaMu2, weight);
      histsIncPresel.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncPresel.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncPresel.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncPresel.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncPresel.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncPresel.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncPresel.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncPresel.nPU->Fill(nPU, weight);
      histsIncPresel.nVtx->Fill(mva.nVtx, weight);
      histsIncPresel.met->Fill(met.pt, weight);

      histsIncPresel.mDiJet->Fill(mva.mDiJet, weight);
      histsIncPresel.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncPresel.ptJet1->Fill(mva.ptJet1, weight);
      histsIncPresel.ptJet2->Fill(mva.ptJet2, weight);
      histsIncPresel.etaJet1->Fill(mva.etaJet1, weight);
      histsIncPresel.etaJet2->Fill(mva.etaJet2, weight);
      histsIncPresel.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncPresel.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncPresel.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncPresel.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncPresel.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncPresel.nJets->Fill(mva.nJets, weight);
      histsIncPresel.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncPresel.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsIncPresel.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsIncPresel.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsIncPresel.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsIncPresel.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsIncPresel.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsIncPresel.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsIncPresel.mDiMuNotBB->Fill(mva.mDiMu, weight);

      if(!vbfPreselection)
      {
        histsIncPresel.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsIncPresel.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsIncPresel.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsIncPresel.BDTHistVBF->Fill(bdtValVBF, weight);
        histsIncPresel.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsIncPresel.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    //Not In Blinding Window Plots
    if (!inBlindWindow)
    {
      histsNotBlindWindow.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsNotBlindWindow.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsNotBlindWindow.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsNotBlindWindow.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsNotBlindWindow.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsNotBlindWindow.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsNotBlindWindow.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsNotBlindWindow.mDiMuNotBB->Fill(mva.mDiMu, weight);

      histsNotBlindWindow.yDiMu->Fill(mva.yDiMu, weight);
      histsNotBlindWindow.ptDiMu->Fill(mva.ptDiMu, weight);
      histsNotBlindWindow.ptMu1->Fill(mva.ptMu1, weight);
      histsNotBlindWindow.ptMu2->Fill(mva.ptMu2, weight);
      histsNotBlindWindow.etaMu1->Fill(mva.etaMu1, weight);
      histsNotBlindWindow.etaMu2->Fill(mva.etaMu2, weight);
      histsNotBlindWindow.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsNotBlindWindow.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsNotBlindWindow.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsNotBlindWindow.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsNotBlindWindow.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsNotBlindWindow.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsNotBlindWindow.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsNotBlindWindow.nPU->Fill(nPU, weight);
      histsNotBlindWindow.nVtx->Fill(mva.nVtx, weight);
      histsNotBlindWindow.met->Fill(met.pt, weight);

      histsNotBlindWindow.mDiJet->Fill(mva.mDiJet, weight);
      histsNotBlindWindow.ptDiJet->Fill(mva.ptDiJet, weight);
      histsNotBlindWindow.ptJet1->Fill(mva.ptJet1, weight);
      histsNotBlindWindow.ptJet2->Fill(mva.ptJet2, weight);
      histsNotBlindWindow.etaJet1->Fill(mva.etaJet1, weight);
      histsNotBlindWindow.etaJet2->Fill(mva.etaJet2, weight);
      histsNotBlindWindow.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsNotBlindWindow.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsNotBlindWindow.deltaRJets->Fill(mva.deltaRJets, weight);
      histsNotBlindWindow.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsNotBlindWindow.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsNotBlindWindow.nJets->Fill(mva.nJets, weight);
      histsNotBlindWindow.ht->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histsNotBlindWindow.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsNotBlindWindow.likelihoodHistMuonOnly->Fill(likeValInc, weight);
      }
      else
      {
        histsNotBlindWindow.BDTHistVBF->Fill(bdtValVBF, weight);
        histsNotBlindWindow.likelihoodHistVBF->Fill(likeValVBF, weight);
      }
    }

    //Upper Control Region Plots
    if (mva.mDiMu>maxBlind)
    {
      histsUpperControlRegion.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsUpperControlRegion.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsUpperControlRegion.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsUpperControlRegion.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsUpperControlRegion.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsUpperControlRegion.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsUpperControlRegion.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsUpperControlRegion.mDiMuNotBB->Fill(mva.mDiMu, weight);

      histsUpperControlRegion.yDiMu->Fill(mva.yDiMu, weight);
      histsUpperControlRegion.ptDiMu->Fill(mva.ptDiMu, weight);
      histsUpperControlRegion.ptMu1->Fill(mva.ptMu1, weight);
      histsUpperControlRegion.ptMu2->Fill(mva.ptMu2, weight);
      histsUpperControlRegion.etaMu1->Fill(mva.etaMu1, weight);
      histsUpperControlRegion.etaMu2->Fill(mva.etaMu2, weight);
      histsUpperControlRegion.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsUpperControlRegion.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsUpperControlRegion.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsUpperControlRegion.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsUpperControlRegion.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsUpperControlRegion.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsUpperControlRegion.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsUpperControlRegion.nPU->Fill(nPU, weight);
      histsUpperControlRegion.nVtx->Fill(mva.nVtx, weight);
      histsUpperControlRegion.met->Fill(met.pt, weight);

      histsUpperControlRegion.mDiJet->Fill(mva.mDiJet, weight);
      histsUpperControlRegion.ptDiJet->Fill(mva.ptDiJet, weight);
      histsUpperControlRegion.ptJet1->Fill(mva.ptJet1, weight);
      histsUpperControlRegion.ptJet2->Fill(mva.ptJet2, weight);
      histsUpperControlRegion.etaJet1->Fill(mva.etaJet1, weight);
      histsUpperControlRegion.etaJet2->Fill(mva.etaJet2, weight);
      histsUpperControlRegion.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsUpperControlRegion.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsUpperControlRegion.deltaRJets->Fill(mva.deltaRJets, weight);
      histsUpperControlRegion.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsUpperControlRegion.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsUpperControlRegion.nJets->Fill(mva.nJets, weight);
      histsUpperControlRegion.ht->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histsUpperControlRegion.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsUpperControlRegion.likelihoodHistMuonOnly->Fill(likeValInc, weight);
      }
      else
      {
        histsUpperControlRegion.BDTHistVBF->Fill(bdtValVBF, weight);
        histsUpperControlRegion.likelihoodHistVBF->Fill(likeValVBF, weight);
      }
    }
    //Lower Control Region Plots
    if (mva.mDiMu<minBlind)
    {
      histsLowerControlRegion.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsLowerControlRegion.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsLowerControlRegion.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsLowerControlRegion.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsLowerControlRegion.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsLowerControlRegion.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsLowerControlRegion.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsLowerControlRegion.mDiMuNotBB->Fill(mva.mDiMu, weight);

      histsLowerControlRegion.yDiMu->Fill(mva.yDiMu, weight);
      histsLowerControlRegion.ptDiMu->Fill(mva.ptDiMu, weight);
      histsLowerControlRegion.ptMu1->Fill(mva.ptMu1, weight);
      histsLowerControlRegion.ptMu2->Fill(mva.ptMu2, weight);
      histsLowerControlRegion.etaMu1->Fill(mva.etaMu1, weight);
      histsLowerControlRegion.etaMu2->Fill(mva.etaMu2, weight);
      histsLowerControlRegion.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsLowerControlRegion.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsLowerControlRegion.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsLowerControlRegion.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsLowerControlRegion.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsLowerControlRegion.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsLowerControlRegion.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsLowerControlRegion.nPU->Fill(nPU, weight);
      histsLowerControlRegion.nVtx->Fill(mva.nVtx, weight);
      histsLowerControlRegion.met->Fill(met.pt, weight);

      histsLowerControlRegion.mDiJet->Fill(mva.mDiJet, weight);
      histsLowerControlRegion.ptDiJet->Fill(mva.ptDiJet, weight);
      histsLowerControlRegion.ptJet1->Fill(mva.ptJet1, weight);
      histsLowerControlRegion.ptJet2->Fill(mva.ptJet2, weight);
      histsLowerControlRegion.etaJet1->Fill(mva.etaJet1, weight);
      histsLowerControlRegion.etaJet2->Fill(mva.etaJet2, weight);
      histsLowerControlRegion.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsLowerControlRegion.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsLowerControlRegion.deltaRJets->Fill(mva.deltaRJets, weight);
      histsLowerControlRegion.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsLowerControlRegion.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsLowerControlRegion.nJets->Fill(mva.nJets, weight);
      histsLowerControlRegion.ht->Fill(mva.ht, weight);
      if(!vbfPreselection)
      {
        histsLowerControlRegion.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsLowerControlRegion.likelihoodHistMuonOnly->Fill(likeValInc, weight);
      }
      else
      {
        histsLowerControlRegion.BDTHistVBF->Fill(bdtValVBF, weight);
        histsLowerControlRegion.likelihoodHistVBF->Fill(likeValVBF, weight);
      }
    }

    if (vbfLoose)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFLoose.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsVBFLoose.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsVBFLoose.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsVBFLoose.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsVBFLoose.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsVBFLoose.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsVBFLoose.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsVBFLoose.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsVBFLoose.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFLoose.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFLoose.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFLoose.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFLoose.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFLoose.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFLoose.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFLoose.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFLoose.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFLoose.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFLoose.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFLoose.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFLoose.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFLoose.nPU->Fill(nPU, weight);
      histsVBFLoose.nVtx->Fill(mva.nVtx, weight);
      histsVBFLoose.met->Fill(met.pt, weight);

      histsVBFLoose.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFLoose.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFLoose.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFLoose.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFLoose.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFLoose.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFLoose.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFLoose.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFLoose.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFLoose.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFLoose.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFLoose.nJets->Fill(mva.nJets, weight);
      histsVBFLoose.ht->Fill(mva.ht, weight);
      histsVBFLoose.BDTHistVBF->Fill(bdtValVBF, weight);
      histsVBFLoose.likelihoodHistVBF->Fill(likeValVBF, weight);
    }

    if (vbfMedium)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFMedium.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsVBFMedium.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsVBFMedium.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsVBFMedium.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsVBFMedium.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsVBFMedium.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsVBFMedium.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsVBFMedium.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsVBFMedium.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFMedium.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFMedium.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFMedium.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFMedium.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFMedium.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFMedium.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFMedium.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFMedium.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFMedium.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFMedium.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFMedium.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFMedium.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFMedium.nPU->Fill(nPU, weight);
      histsVBFMedium.nVtx->Fill(mva.nVtx, weight);
      histsVBFMedium.met->Fill(met.pt, weight);

      histsVBFMedium.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFMedium.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFMedium.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFMedium.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFMedium.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFMedium.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFMedium.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFMedium.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFMedium.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFMedium.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFMedium.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFMedium.nJets->Fill(mva.nJets, weight);
      histsVBFMedium.ht->Fill(mva.ht, weight);
      histsVBFMedium.BDTHistVBF->Fill(bdtValVBF, weight);
      histsVBFMedium.likelihoodHistVBF->Fill(likeValVBF, weight);
    }

    if (vbfTight)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFTight.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsVBFTight.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsVBFTight.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsVBFTight.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsVBFTight.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsVBFTight.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsVBFTight.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsVBFTight.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsVBFTight.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFTight.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFTight.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFTight.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFTight.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFTight.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFTight.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFTight.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFTight.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFTight.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFTight.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFTight.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFTight.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFTight.nPU->Fill(nPU, weight);
      histsVBFTight.nVtx->Fill(mva.nVtx, weight);
      histsVBFTight.met->Fill(met.pt, weight);

      histsVBFTight.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFTight.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFTight.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFTight.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFTight.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFTight.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFTight.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFTight.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFTight.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFTight.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFTight.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFTight.nJets->Fill(mva.nJets, weight);
      histsVBFTight.ht->Fill(mva.ht, weight);
      histsVBFTight.BDTHistVBF->Fill(bdtValVBF, weight);
      histsVBFTight.likelihoodHistVBF->Fill(likeValVBF, weight);
    }

    if (vbfVeryTight)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFVeryTight.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsVBFVeryTight.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsVBFVeryTight.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsVBFVeryTight.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsVBFVeryTight.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsVBFVeryTight.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsVBFVeryTight.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsVBFVeryTight.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsVBFVeryTight.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFVeryTight.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFVeryTight.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFVeryTight.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFVeryTight.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFVeryTight.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFVeryTight.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFVeryTight.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFVeryTight.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFVeryTight.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFVeryTight.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFVeryTight.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFVeryTight.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFVeryTight.nPU->Fill(nPU, weight);
      histsVBFVeryTight.nVtx->Fill(mva.nVtx, weight);
      histsVBFVeryTight.met->Fill(met.pt, weight);

      histsVBFVeryTight.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFVeryTight.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFVeryTight.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFVeryTight.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFVeryTight.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFVeryTight.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFVeryTight.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFVeryTight.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFVeryTight.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFVeryTight.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFVeryTight.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFVeryTight.nJets->Fill(mva.nJets, weight);
      histsVBFVeryTight.ht->Fill(mva.ht, weight);
      histsVBFVeryTight.BDTHistVBF->Fill(bdtValVBF, weight);
      histsVBFVeryTight.likelihoodHistVBF->Fill(likeValVBF, weight);
    }

    if (pt0to30)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsPt0to30.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsPt0to30.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsPt0to30.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsPt0to30.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsPt0to30.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsPt0to30.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsPt0to30.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsPt0to30.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsPt0to30.yDiMu->Fill(mva.yDiMu, weight);
      histsPt0to30.ptDiMu->Fill(mva.ptDiMu, weight);
      histsPt0to30.ptMu1->Fill(mva.ptMu1, weight);
      histsPt0to30.ptMu2->Fill(mva.ptMu2, weight);
      histsPt0to30.etaMu1->Fill(mva.etaMu1, weight);
      histsPt0to30.etaMu2->Fill(mva.etaMu2, weight);
      histsPt0to30.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsPt0to30.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsPt0to30.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsPt0to30.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsPt0to30.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsPt0to30.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsPt0to30.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsPt0to30.nPU->Fill(nPU, weight);
      histsPt0to30.nVtx->Fill(mva.nVtx, weight);
      histsPt0to30.met->Fill(met.pt, weight);

      histsPt0to30.mDiJet->Fill(mva.mDiJet, weight);
      histsPt0to30.ptDiJet->Fill(mva.ptDiJet, weight);
      histsPt0to30.ptJet1->Fill(mva.ptJet1, weight);
      histsPt0to30.ptJet2->Fill(mva.ptJet2, weight);
      histsPt0to30.etaJet1->Fill(mva.etaJet1, weight);
      histsPt0to30.etaJet2->Fill(mva.etaJet2, weight);
      histsPt0to30.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsPt0to30.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsPt0to30.deltaRJets->Fill(mva.deltaRJets, weight);
      histsPt0to30.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsPt0to30.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsPt0to30.nJets->Fill(mva.nJets, weight);
      histsPt0to30.ht->Fill(mva.ht, weight);
      histsPt0to30.BDTHistMuonOnly->Fill(bdtValInc, weight);
      histsPt0to30.likelihoodHistMuonOnly->Fill(likeValInc, weight);
    }

    if (pt30to50)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsPt30to50.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsPt30to50.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsPt30to50.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsPt30to50.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsPt30to50.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsPt30to50.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsPt30to50.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsPt30to50.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsPt30to50.yDiMu->Fill(mva.yDiMu, weight);
      histsPt30to50.ptDiMu->Fill(mva.ptDiMu, weight);
      histsPt30to50.ptMu1->Fill(mva.ptMu1, weight);
      histsPt30to50.ptMu2->Fill(mva.ptMu2, weight);
      histsPt30to50.etaMu1->Fill(mva.etaMu1, weight);
      histsPt30to50.etaMu2->Fill(mva.etaMu2, weight);
      histsPt30to50.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsPt30to50.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsPt30to50.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsPt30to50.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsPt30to50.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsPt30to50.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsPt30to50.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsPt30to50.nPU->Fill(nPU, weight);
      histsPt30to50.nVtx->Fill(mva.nVtx, weight);
      histsPt30to50.met->Fill(met.pt, weight);

      histsPt30to50.mDiJet->Fill(mva.mDiJet, weight);
      histsPt30to50.ptDiJet->Fill(mva.ptDiJet, weight);
      histsPt30to50.ptJet1->Fill(mva.ptJet1, weight);
      histsPt30to50.ptJet2->Fill(mva.ptJet2, weight);
      histsPt30to50.etaJet1->Fill(mva.etaJet1, weight);
      histsPt30to50.etaJet2->Fill(mva.etaJet2, weight);
      histsPt30to50.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsPt30to50.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsPt30to50.deltaRJets->Fill(mva.deltaRJets, weight);
      histsPt30to50.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsPt30to50.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsPt30to50.nJets->Fill(mva.nJets, weight);
      histsPt30to50.ht->Fill(mva.ht, weight);
      histsPt30to50.BDTHistMuonOnly->Fill(bdtValInc, weight);
      histsPt30to50.likelihoodHistMuonOnly->Fill(likeValInc, weight);
    }

    if (pt50to125)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsPt50to125.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsPt50to125.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsPt50to125.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsPt50to125.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsPt50to125.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsPt50to125.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsPt50to125.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsPt50to125.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsPt50to125.yDiMu->Fill(mva.yDiMu, weight);
      histsPt50to125.ptDiMu->Fill(mva.ptDiMu, weight);
      histsPt50to125.ptMu1->Fill(mva.ptMu1, weight);
      histsPt50to125.ptMu2->Fill(mva.ptMu2, weight);
      histsPt50to125.etaMu1->Fill(mva.etaMu1, weight);
      histsPt50to125.etaMu2->Fill(mva.etaMu2, weight);
      histsPt50to125.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsPt50to125.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsPt50to125.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsPt50to125.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsPt50to125.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsPt50to125.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsPt50to125.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsPt50to125.nPU->Fill(nPU, weight);
      histsPt50to125.nVtx->Fill(mva.nVtx, weight);
      histsPt50to125.met->Fill(met.pt, weight);

      histsPt50to125.mDiJet->Fill(mva.mDiJet, weight);
      histsPt50to125.ptDiJet->Fill(mva.ptDiJet, weight);
      histsPt50to125.ptJet1->Fill(mva.ptJet1, weight);
      histsPt50to125.ptJet2->Fill(mva.ptJet2, weight);
      histsPt50to125.etaJet1->Fill(mva.etaJet1, weight);
      histsPt50to125.etaJet2->Fill(mva.etaJet2, weight);
      histsPt50to125.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsPt50to125.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsPt50to125.deltaRJets->Fill(mva.deltaRJets, weight);
      histsPt50to125.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsPt50to125.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsPt50to125.nJets->Fill(mva.nJets, weight);
      histsPt50to125.ht->Fill(mva.ht, weight);
      histsPt50to125.BDTHistMuonOnly->Fill(bdtValInc, weight);
      histsPt50to125.likelihoodHistMuonOnly->Fill(likeValInc, weight);
    }

    if (pt125to250)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsPt125to250.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsPt125to250.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsPt125to250.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsPt125to250.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsPt125to250.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsPt125to250.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsPt125to250.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsPt125to250.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsPt125to250.yDiMu->Fill(mva.yDiMu, weight);
      histsPt125to250.ptDiMu->Fill(mva.ptDiMu, weight);
      histsPt125to250.ptMu1->Fill(mva.ptMu1, weight);
      histsPt125to250.ptMu2->Fill(mva.ptMu2, weight);
      histsPt125to250.etaMu1->Fill(mva.etaMu1, weight);
      histsPt125to250.etaMu2->Fill(mva.etaMu2, weight);
      histsPt125to250.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsPt125to250.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsPt125to250.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsPt125to250.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsPt125to250.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsPt125to250.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsPt125to250.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsPt125to250.nPU->Fill(nPU, weight);
      histsPt125to250.nVtx->Fill(mva.nVtx, weight);
      histsPt125to250.met->Fill(met.pt, weight);

      histsPt125to250.mDiJet->Fill(mva.mDiJet, weight);
      histsPt125to250.ptDiJet->Fill(mva.ptDiJet, weight);
      histsPt125to250.ptJet1->Fill(mva.ptJet1, weight);
      histsPt125to250.ptJet2->Fill(mva.ptJet2, weight);
      histsPt125to250.etaJet1->Fill(mva.etaJet1, weight);
      histsPt125to250.etaJet2->Fill(mva.etaJet2, weight);
      histsPt125to250.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsPt125to250.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsPt125to250.deltaRJets->Fill(mva.deltaRJets, weight);
      histsPt125to250.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsPt125to250.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsPt125to250.nJets->Fill(mva.nJets, weight);
      histsPt125to250.ht->Fill(mva.ht, weight);
      histsPt125to250.BDTHistMuonOnly->Fill(bdtValInc, weight);
      histsPt125to250.likelihoodHistMuonOnly->Fill(likeValInc, weight);
    }

    if (pt250)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsPt250.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsPt250.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsPt250.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsPt250.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsPt250.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsPt250.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsPt250.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsPt250.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsPt250.yDiMu->Fill(mva.yDiMu, weight);
      histsPt250.ptDiMu->Fill(mva.ptDiMu, weight);
      histsPt250.ptMu1->Fill(mva.ptMu1, weight);
      histsPt250.ptMu2->Fill(mva.ptMu2, weight);
      histsPt250.etaMu1->Fill(mva.etaMu1, weight);
      histsPt250.etaMu2->Fill(mva.etaMu2, weight);
      histsPt250.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsPt250.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsPt250.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsPt250.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsPt250.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsPt250.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsPt250.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsPt250.nPU->Fill(nPU, weight);
      histsPt250.nVtx->Fill(mva.nVtx, weight);
      histsPt250.met->Fill(met.pt, weight);

      histsPt250.mDiJet->Fill(mva.mDiJet, weight);
      histsPt250.ptDiJet->Fill(mva.ptDiJet, weight);
      histsPt250.ptJet1->Fill(mva.ptJet1, weight);
      histsPt250.ptJet2->Fill(mva.ptJet2, weight);
      histsPt250.etaJet1->Fill(mva.etaJet1, weight);
      histsPt250.etaJet2->Fill(mva.etaJet2, weight);
      histsPt250.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsPt250.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsPt250.deltaRJets->Fill(mva.deltaRJets, weight);
      histsPt250.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsPt250.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsPt250.nJets->Fill(mva.nJets, weight);
      histsPt250.ht->Fill(mva.ht, weight);
      histsPt250.BDTHistMuonOnly->Fill(bdtValInc, weight);
      histsPt250.likelihoodHistMuonOnly->Fill(likeValInc, weight);
    }

    if (!vbfPreselection && mva.getMVAPassBDTCut(cfgNameInc))
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncBDTSig80.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsIncBDTSig80.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsIncBDTSig80.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsIncBDTSig80.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsIncBDTSig80.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsIncBDTSig80.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsIncBDTSig80.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsIncBDTSig80.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsIncBDTSig80.yDiMu->Fill(mva.yDiMu, weight);
      histsIncBDTSig80.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncBDTSig80.ptMu1->Fill(mva.ptMu1, weight);
      histsIncBDTSig80.ptMu2->Fill(mva.ptMu2, weight);
      histsIncBDTSig80.etaMu1->Fill(mva.etaMu1, weight);
      histsIncBDTSig80.etaMu2->Fill(mva.etaMu2, weight);
      histsIncBDTSig80.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncBDTSig80.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncBDTSig80.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncBDTSig80.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncBDTSig80.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncBDTSig80.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncBDTSig80.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncBDTSig80.nPU->Fill(nPU, weight);
      histsIncBDTSig80.nVtx->Fill(mva.nVtx, weight);
      histsIncBDTSig80.met->Fill(met.pt, weight);

      histsIncBDTSig80.mDiJet->Fill(mva.mDiJet, weight);
      histsIncBDTSig80.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncBDTSig80.ptJet1->Fill(mva.ptJet1, weight);
      histsIncBDTSig80.ptJet2->Fill(mva.ptJet2, weight);
      histsIncBDTSig80.etaJet1->Fill(mva.etaJet1, weight);
      histsIncBDTSig80.etaJet2->Fill(mva.etaJet2, weight);
      histsIncBDTSig80.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncBDTSig80.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncBDTSig80.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncBDTSig80.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncBDTSig80.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncBDTSig80.nJets->Fill(mva.nJets, weight);
      histsIncBDTSig80.ht->Fill(mva.ht, weight);
      histsIncBDTSig80.BDTHistMuonOnly->Fill(bdtValInc, weight);
    }

    if (vbfPreselection && mva.getMVAPassBDTCut(cfgNameVBF))
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFBDTSig80.mDiMu->Fill(mva.mDiMu, weight);
      if (isBB) histsVBFBDTSig80.mDiMuBB->Fill(mva.mDiMu, weight);
      if (isBE) histsVBFBDTSig80.mDiMuBE->Fill(mva.mDiMu, weight);
      if (isBO) histsVBFBDTSig80.mDiMuBO->Fill(mva.mDiMu, weight);
      if (isOO) histsVBFBDTSig80.mDiMuOO->Fill(mva.mDiMu, weight);
      if (isOE) histsVBFBDTSig80.mDiMuOE->Fill(mva.mDiMu, weight);
      if (isEE) histsVBFBDTSig80.mDiMuEE->Fill(mva.mDiMu, weight);
      if (isNotBB) histsVBFBDTSig80.mDiMuNotBB->Fill(mva.mDiMu, weight);

#ifdef BLIND
      }
#endif
      histsVBFBDTSig80.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFBDTSig80.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFBDTSig80.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFBDTSig80.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFBDTSig80.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFBDTSig80.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFBDTSig80.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFBDTSig80.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFBDTSig80.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFBDTSig80.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFBDTSig80.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFBDTSig80.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFBDTSig80.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFBDTSig80.nPU->Fill(nPU, weight);
      histsVBFBDTSig80.nVtx->Fill(mva.nVtx, weight);
      histsVBFBDTSig80.met->Fill(met.pt, weight);

      histsVBFBDTSig80.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFBDTSig80.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFBDTSig80.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFBDTSig80.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFBDTSig80.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFBDTSig80.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFBDTSig80.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFBDTSig80.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFBDTSig80.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFBDTSig80.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFBDTSig80.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFBDTSig80.nJets->Fill(mva.nJets, weight);
      histsVBFBDTSig80.ht->Fill(mva.ht, weight);
      histsVBFBDTSig80.BDTHistVBF->Fill(bdtValVBF, weight);
    }

    timeReading += difftime(timeStopReading,timeStartReading);
    timeProcessing += difftime(timeStartFilling,timeStopReading);
    timeFilling += difftime(time(NULL),timeStartFilling);
  }// end event loop
  time_t timeEndEventLoop = time(NULL);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  outFile->cd();

  hists.Write();

  TDirectory* dirCalibUp = outFile->mkdir("CalibUp");
  dirCalibUp->cd();
  histsCalibUp.Write();

  TDirectory* dirCalibDown = outFile->mkdir("CalibDown");
  dirCalibDown->cd();
  histsCalibDown.Write();

  TDirectory* dirResUp = outFile->mkdir("ResUp");
  dirResUp->cd();
  histsResUp.Write();

  TDirectory* dirResDown = outFile->mkdir("ResDown");
  dirResDown->cd();
  histsResDown.Write();

  TDirectory* dir4GeVWindow = outFile->mkdir("4GeVWindow");
  dir4GeVWindow->cd();
  hists4GeVWindow.Write();

  TDirectory* dirPtDiMu100 = outFile->mkdir("PtDiMu100");
  dirPtDiMu100->cd();
  histsPtDiMu100.Write();

  TDirectory* dirVBFPresel = outFile->mkdir("VBFPresel");
  dirVBFPresel->cd();
  histsVBFPresel.Write();

  TDirectory* dirIncPresel = outFile->mkdir("IncPresel");
  dirIncPresel->cd();
  histsIncPresel.Write();

  TDirectory* dirNotBlindWindow = outFile->mkdir("NotBlindWindow");
  dirNotBlindWindow->cd();
  histsNotBlindWindow.Write();

  TDirectory* dirUpperControlRegion = outFile->mkdir("UpperControlRegion");
  dirUpperControlRegion->cd();
  histsUpperControlRegion.Write();

  TDirectory* dirLowerControlRegion = outFile->mkdir("LowerControlRegion");
  dirLowerControlRegion->cd();
  histsLowerControlRegion.Write();

  TDirectory* dirVBFLoose = outFile->mkdir("VBFLoose");
  dirVBFLoose->cd();
  histsVBFLoose.Write();

  TDirectory* dirVBFMedium = outFile->mkdir("VBFMedium");
  dirVBFMedium->cd();
  histsVBFMedium.Write();

  TDirectory* dirVBFTight = outFile->mkdir("VBFTight");
  dirVBFTight->cd();
  histsVBFTight.Write();

  TDirectory* dirVBFVeryTight = outFile->mkdir("VBFVeryTight");
  dirVBFVeryTight->cd();
  histsVBFVeryTight.Write();

  TDirectory* dirPt0to30 = outFile->mkdir("Pt0to30");
  dirPt0to30->cd();
  histsPt0to30.Write();

  TDirectory* dirPt30to50 = outFile->mkdir("Pt30to50");
  dirPt30to50->cd();
  histsPt30to50.Write();

  TDirectory* dirPt50to125 = outFile->mkdir("Pt50to125");
  dirPt50to125->cd();
  histsPt50to125.Write();

  TDirectory* dirPt125to250 = outFile->mkdir("Pt125to250");
  dirPt125to250->cd();
  histsPt125to250.Write();

  TDirectory* dirPt250 = outFile->mkdir("Pt250");
  dirPt250->cd();
  histsPt250.Write();

  TDirectory* dirIncBDTSig80 = outFile->mkdir("IncBDTSig80");
  dirIncBDTSig80->cd();
  histsIncBDTSig80.Write();

  TDirectory* dirVBFBDTSig80 = outFile->mkdir("VBFBDTSig80");
  dirVBFBDTSig80->cd();
  histsVBFBDTSig80.Write();

  ofstream testOutFile;
  testOutFile.open("testEventNums.txt");
  testOutFile << testString;
  testOutFile.close();
  cout <<"#######################\n"<< testString <<"#######################\n" << endl;
  cout << "testCounter: "<< testCounter << endl;
  cout << "Total Time: "<<std::setprecision(3) << difftime(time(NULL),timeStart)<<"\n";
  cout << "Setup Time: "<<std::setprecision(3) <<difftime(timeStartEventLoop,timeStart)<<"\n";
  cout << "Event Loop Time: "<<std::setprecision(3) 
        <<difftime(timeEndEventLoop,timeStartEventLoop)<< ", "<<std::setprecision(3) 
        <<difftime(timeEndEventLoop,timeStartEventLoop)/(std::min(nEvents,(unsigned) maxEvents))*1000.
        <<" s / 1000 events or "
        <<(std::min(nEvents,(unsigned) maxEvents))/difftime(timeEndEventLoop,timeStartEventLoop)*3600./1.0e6
        <<"M events/hour \n";
  cout << "  Read Time: "<<std::setprecision(3) << timeReading << std::endl;
  cout << "  Proc Time: "<<std::setprecision(3) << timeProcessing << std::endl;
  cout << "  Fill Time: "<<std::setprecision(3) << timeFilling << std::endl;
  cout << "  All Read Time: "<<std::setprecision(3) << timeReadingAll << std::endl;
  cout << "Wrapup Time: "<<std::setprecision(3) <<difftime(time(NULL),timeEndEventLoop)<<"\n";
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

HistStruct::HistStruct()
{

  unsigned nMassBins = 400;
  float minMass = 0.;
  float maxMass = 400.;
  unsigned nMVABins = 200;
  mDiMu = new TH1F("mDiMu","DiMuon Mass",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMu);

  mDiMuBB = new TH1F("mDiMuBB","DiMuon Mass Barrel-Barrel",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuBB);
  mDiMuBO = new TH1F("mDiMuBO","DiMuon Mass Barrel-Overlap",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuBO);
  mDiMuBE = new TH1F("mDiMuBE","DiMuon Mass Barrel-Endcap",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuBE);
  mDiMuOO = new TH1F("mDiMuOO","DiMuon Mass Overlap-Overlap",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuOO);
  mDiMuOE = new TH1F("mDiMuOE","DiMuon Mass Overlap-Endcap",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuOE);
  mDiMuEE = new TH1F("mDiMuEE","DiMuon Mass Endcap-Endcap",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuEE);
  mDiMuNotBB = new TH1F("mDiMuNotBB","DiMuon Mass Not Barrel-Barrel",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuNotBB);

  mDiJet = new TH1F("mDiJet","DiJet Mass",500,0,2000);
  histVec.push_back(mDiJet);

  ptDiMu = new TH1F("ptDiMu","DiMuon Pt",250,0,500);
  histVec.push_back(ptDiMu);

  ptDiJet = new TH1F("ptDiJet","DiJet Pt",250,0,1000);
  histVec.push_back(ptDiJet);

  yDiMu = new TH1F("yDiMu","DiMuon Rapidity",100,-4,4);
  histVec.push_back(yDiMu);

  yVptDiMu = new TH2F("yVptDiMu","DiMuon Rapidity v. p_{T}",250,0,500,100,0,4);
  histVec2D.push_back(yVptDiMu);
  ptVmDiMu = new TH2F("ptVmDiMu","DiMuon p_{T} v. Mass",nMassBins,minMass,maxMass,250,0,250);
  histVec2D.push_back(ptVmDiMu);
  yVmDiMu = new TH2F("yVmDiMu","DiMuon |y| v. Mass",nMassBins,minMass,maxMass,100,0,4);
  phiVmDiMu = new TH2F("phiVmDiMu","DiMuon #phi v. Mass",nMassBins,minMass,maxMass,100,0,3.2);
  histVec2D.push_back(phiVmDiMu);

  ptMu1 = new TH1F("ptMu1","Leading Muon Pt",100,0,200);
  histVec.push_back(ptMu1);
  ptMu2 = new TH1F("ptMu2","Sub-Leading Muon Pt",100,0,200);
  histVec.push_back(ptMu2);
  ptJet1 = new TH1F("ptJet1","Leading Jet Pt",200,0,1000);
  histVec.push_back(ptJet1);
  ptJet2 = new TH1F("ptJet2","Sub-Leading Jet Pt",200,0,1000);
  histVec.push_back(ptJet2);

  etaMu1 = new TH1F("etaMu1","Leading Muon #eta",50,-2.5,2.5);
  histVec.push_back(etaMu1);
  etaMu2 = new TH1F("etaMu2","Sub-Leading Muon #eta",50,-2.5,2.5);
  histVec.push_back(etaMu2);
  etaJet1 = new TH1F("etaJet1","Leading Jet #eta",50,-5.0,5.0);
  histVec.push_back(etaJet1);
  etaJet2 = new TH1F("etaJet2","Sub-Leading Jet #eta",50,-5.0,5.0);
  histVec.push_back(etaJet2);

  deltaEtaJets = new TH1F("deltaEtaJets","#Delta#eta Jets",50,0.0,10.0);
  histVec.push_back(deltaEtaJets);
  deltaPhiJets = new TH1F("deltaPhiJets","#Delta#phi Jets",50,0.0,3.2);
  histVec.push_back(deltaPhiJets);
  deltaRJets = new TH1F("deltaRJets","#Delta R Jets",50,0.0,10.0);
  histVec.push_back(deltaRJets);

  deltaEtaMuons = new TH1F("deltaEtaMuons","#Delta#eta Jets",50,0.0,10.0);
  histVec.push_back(deltaEtaMuons);
  deltaPhiMuons = new TH1F("deltaPhiMuons","#Delta#phi Jets",50,0.0,3.2);
  histVec.push_back(deltaPhiMuons);
  deltaRMuons = new TH1F("deltaRMuons","#Delta R Jets",50,0.0,10.0);
  histVec.push_back(deltaRMuons);

  countsHist = new TH1F("countsHist","Event Counts",10,0.0,10.0);
  countsHist->GetXaxis()->SetBinLabel(1,"total");
  countsHist->GetXaxis()->SetBinLabel(2,"2#mu ID");
  countsHist->GetXaxis()->SetBinLabel(3,"HLT");
  countsHist->GetXaxis()->SetBinLabel(4,"Charge");
  countsHist->GetXaxis()->SetBinLabel(5,"m_{#mu#mu}");
  countsHist->GetXaxis()->SetBinLabel(6,"Inc Pre");
  countsHist->GetXaxis()->SetBinLabel(7,"VBF Pre");
  histVec.push_back(countsHist);

  countsHist2 = new TH1F("countsHist2","Event Counts",14,0.0,14.0);
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
  histVec.push_back(countsHist2);

  cosThetaStar = new TH1F("cosThetaStar","cos(#theta^{*})",50,-1.,1.);
  histVec.push_back(cosThetaStar);
  cosThetaStarCS = new TH1F("cosThetaStarCS","cos(#theta^{*}_{CS})",50,-1.,1.);
  histVec.push_back(cosThetaStarCS);

  puJetIDSimpleDiscJet1 = new TH1F("puJetIDSimpleDiscJet1","PU Jet ID--Simple Discriminator Leading Jet",50,-1.,1.);
  histVec.push_back(puJetIDSimpleDiscJet1);
  puJetIDSimpleDiscJet2 = new TH1F("puJetIDSimpleDiscJet2","PU Jet ID--Simple Discriminator Sub-Leading Jet",50,-1.,1.);
  histVec.push_back(puJetIDSimpleDiscJet2);
  puJetIDSimpleDiscJet3 = new TH1F("puJetIDSimpleDiscJet3","PU Jet ID--Simple Discriminator 3rd Leading Jet",50,-1.,1.);
  histVec.push_back(puJetIDSimpleDiscJet3);

  puJetIDSimpleJet1 = new TH1F("puJetIDSimpleJet1","PU Jet ID--Simple Loose Leading Jet",2,0,2);
  histVec.push_back(puJetIDSimpleJet1);
  puJetIDSimpleJet1->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet1->GetXaxis()->SetBinLabel(2,"Pass");
  puJetIDSimpleJet2 = new TH1F("puJetIDSimpleJet2","PU Jet ID--Simple Loose Sub-Leading Jet",2,-0,2);
  histVec.push_back(puJetIDSimpleJet2);
  puJetIDSimpleJet2->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet2->GetXaxis()->SetBinLabel(2,"Pass");
  puJetIDSimpleJet3 = new TH1F("puJetIDSimpleJet3","PU Jet ID--Simple Loose 3rd Leading Jet",2,-0,2);
  histVec.push_back(puJetIDSimpleJet3);
  puJetIDSimpleJet3->GetXaxis()->SetBinLabel(1,"Fail");
  puJetIDSimpleJet3->GetXaxis()->SetBinLabel(2,"Pass");

  BDTHistMuonOnly = new TH1F("BDTHistMuonOnly","BDT Discriminator",nMVABins,-1,1);
  histVec.push_back(BDTHistMuonOnly);
  likelihoodHistMuonOnly = new TH1F("likelihoodHistMuonOnly","Likelihood Discriminator",nMVABins,-1,1);
  histVec.push_back(likelihoodHistMuonOnly);

  BDTHistVBF = new TH1F("BDTHistVBF","BDT Discriminator",nMVABins,-1,1);
  histVec.push_back(BDTHistVBF);
  likelihoodHistVBF = new TH1F("likelihoodHistVBF","Likelihood Discriminator",nMVABins,-1,1);
  histVec.push_back(likelihoodHistVBF);

  BDTHistMuonOnlyVMass = new TH2F("BDTHistMuonOnlyVMass","BDT Discriminator",nMassBins,minMass,maxMass,nMVABins,-1,1);
  histVec2D.push_back(BDTHistMuonOnlyVMass);
  likelihoodHistMuonOnlyVMass = new TH2F("likelihoodHistMuonOnlyVMass","Likelihood Discriminator",nMassBins,minMass,maxMass,nMVABins,-1,1);
  histVec2D.push_back(likelihoodHistMuonOnlyVMass);
  BDTHistVBFVMass = new TH2F("BDTHistVBFVMass","BDT Discriminator",nMassBins,minMass,maxMass,nMVABins,-1,1);
  histVec2D.push_back(BDTHistVBFVMass);
  likelihoodHistVBFVMass = new TH2F("likelihoodHistVBFVMass","Likelihood Discriminator",nMassBins,minMass,maxMass,nMVABins,-1,1);
  histVec2D.push_back(likelihoodHistVBFVMass);

  relIsoMu1 = new TH1F("relIsoMu1","",1000,0,10.0);
  histVec.push_back(relIsoMu1);
  relIsoMu2 = new TH1F("relIsoMu2","",1000,0,10.0);
  histVec.push_back(relIsoMu2);

  nJets = new TH1F("nJets","",11,0,11);
  histVec.push_back(nJets);
  ht = new TH1F("ht","",200,0,2000);
  histVec.push_back(ht);
  nJetsInRapidityGap = new TH1F("nJetsInRapidityGap","",11,0,11);
  histVec.push_back(nJetsInRapidityGap);
  htInRapidityGap = new TH1F("htInRapidityGap","",200,0,2000);
  histVec.push_back(htInRapidityGap);

  nPU = new TH1F("nPU","",100,0,100);
  histVec.push_back(nPU);
  nVtx = new TH1F("nVtx","",100,0,100);
  histVec.push_back(nVtx);
  met = new TH1F("met","",160,0,800);
  histVec.push_back(met);
  weight = new TH1F("weight","",500,0,5.0);
  histVec.push_back(weight);

  std::vector<TH1F*>::iterator hist;
  std::vector<TH2F*>::iterator hist2D;
  for(hist = histVec.begin();hist != histVec.end(); hist++)
    (*hist)->Sumw2();
  for(hist2D = histVec2D.begin();hist2D != histVec2D.end(); hist2D++)
    (*hist2D)->Sumw2();
}

HistStruct::~HistStruct()
{
  std::vector<TH1F*>::iterator hist;
  std::vector<TH2F*>::iterator hist2D;
  for(hist = histVec.begin();hist != histVec.end(); hist++)
    delete *hist;
  for(hist2D = histVec2D.begin();hist2D != histVec2D.end(); hist2D++)
    delete *hist2D;
}

void
HistStruct::Write()
{
  std::vector<TH1F*>::iterator hist;
  std::vector<TH2F*>::iterator hist2D;
  for(hist = histVec.begin();hist != histVec.end(); hist++)
    (*hist)->Write();
  for(hist2D = histVec2D.begin();hist2D != histVec2D.end(); hist2D++)
    (*hist2D)->Write();
}
