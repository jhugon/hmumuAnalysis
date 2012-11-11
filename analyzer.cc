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

#include "annaCalibCode/SmearingTool.h"

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
  void Write(TFile* outfile, std::string directory);

  std::vector<TH1F*> histVec;
  std::vector<TH2F*> histVec2D;

  TH1F* mDiMu;
  TH1F* mDiMuResSigUp;
  TH1F* mDiMuResSigDown;
  TH1F* mDiMuResASigUp;
  TH1F* mDiMuResASigDown;

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
        break;
    }
  }

  // Check to see if it is signal
  bool isSignal = false;
  std::vector<std::string> signalWords;
  signalWords.push_back("Hmumu");
  std::vector<std::string>::const_iterator signalWord;
  for(signalWord = signalWords.begin(); signalWord != signalWords.end();signalWord++)
  {
    regex re(*signalWord);
    if(regex_search(filenames[0],re))
    {
        isSignal = true;
        break;
    }
  }

  if(isData)
    std::cout << "This is a Real Data Sample\n";
  else
    std::cout << "This is a MC Sample\n";
  if(isSignal)
    std::cout << "This is a Signal Sample\n";

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

  _TrackInfo reco1GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM1HpostFSR"))
    tree->SetBranchAddress("genM1HpostFSR", &reco1GenPostFSR);

  _TrackInfo reco2GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM2HpostFSR"))
    tree->SetBranchAddress("genM2HpostFSR", &reco2GenPostFSR);

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
  HistStruct histsBB;
  HistStruct histsBO;
  HistStruct histsBE;
  HistStruct histsOO;
  HistStruct histsOE;
  HistStruct histsEE;
  HistStruct histsNotBB;

  HistStruct hists4GeVWindow;

  HistStruct histsIncPresel;
  HistStruct histsIncPreselBB;
  HistStruct histsIncPreselBO;
  HistStruct histsIncPreselBE;
  HistStruct histsIncPreselOO;
  HistStruct histsIncPreselOE;
  HistStruct histsIncPreselEE;
  HistStruct histsIncPreselNotBB;

  HistStruct histsVBFPresel;
  HistStruct histsVBFPreselBB;
  HistStruct histsVBFPreselNotBB;

  HistStruct histsIncBDTSig80;
  HistStruct histsIncBDTSig80BB;
  HistStruct histsIncBDTSig80BO;
  HistStruct histsIncBDTSig80BE;
  HistStruct histsIncBDTSig80OO;
  HistStruct histsIncBDTSig80OE;
  HistStruct histsIncBDTSig80EE;
  HistStruct histsIncBDTSig80NotBB;

  HistStruct histsVBFBDTSig80;
  HistStruct histsVBFBDTSig80BB;
  HistStruct histsVBFBDTSig80NotBB;

  HistStruct histsIncPreselDiMuPtL20;
  HistStruct histsVBFPreselDiMuPtL20;

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

  /////////////////////////
  // Smearing
  SmearingTool *smearPT = new SmearingTool();

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

    float mDiMuResSigUp = recoCandMass;
    float mDiMuResSigDown = recoCandMass;
    float mDiMuResASigUp = recoCandMass;
    float mDiMuResASigDown = recoCandMass;

#ifdef SMEARING
    if(isSignal)
    {
      if(reco1GenPostFSR.pt<0.)
        cout << "Muon 1 Post FSR not valid!\n";
      if(reco2GenPostFSR.pt<0.)
        cout << "Muon 2 Post FSR not valid!\n";
      float ptReco1 = smearPT -> PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge);
      float ptReco2 = smearPT -> PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge);
      TLorentzVector reco1Vec;
      TLorentzVector reco2Vec;
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      TLorentzVector diMuonVec = reco1Vec + reco2Vec;

      reco1.pt = ptReco1;
      reco2.pt = ptReco2;
      recoCandMass = diMuonVec.M();
      recoCandPt = diMuonVec.Pt();
      recoCandY = diMuonVec.Rapidity();
      recoCandPhi = diMuonVec.Phi();

      //Systematics Time
      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge,"sig1",1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge,"sig1",1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResSigUp = diMuonVec.M();

      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge,"sig1",-1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge,"sig1",-1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResSigDown = diMuonVec.M();

      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge,"Asig2Var",1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge,"Asig2Var",1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResASigUp = diMuonVec.M();

      ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge,"Asig2Var",-1.0);
      ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge,"Asig2Var",-1.0);
      reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,0.105);
      reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,0.105);
      diMuonVec = reco1Vec + reco2Vec;
      mDiMuResASigDown = diMuonVec.M();
      
    }
#endif

    fillMuonHist(hists.countsHist2, reco1, reco2);
    //printStationMiss(reco1,reco2,eventInfo,testString,testCounter);

    mva.resetValues();
    mva.mDiMu = recoCandMass;
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
    hists.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
    hists.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
    hists.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
    hists.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);
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

    }
    else
    {
      hists.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      hists.likelihoodHistVBFVMass->Fill(mva.mDiMu, likeValVBF, weight);
      hists.BDTHistVBF->Fill(bdtValVBF, weight);
      hists.likelihoodHistVBF->Fill(likeValVBF, weight);

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
      hists4GeVWindow.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      hists4GeVWindow.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      hists4GeVWindow.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      hists4GeVWindow.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

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

    if (isBB)
    {
      histsBB.yDiMu->Fill(mva.yDiMu, weight);
      histsBB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsBB.ptMu1->Fill(mva.ptMu1, weight);
      histsBB.ptMu2->Fill(mva.ptMu2, weight);
      histsBB.etaMu1->Fill(mva.etaMu1, weight);
      histsBB.etaMu2->Fill(mva.etaMu2, weight);
      histsBB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsBB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsBB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsBB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsBB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsBB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsBB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsBB.nPU->Fill(nPU, weight);
      histsBB.nVtx->Fill(mva.nVtx, weight);
      histsBB.met->Fill(met.pt, weight);

      histsBB.mDiJet->Fill(mva.mDiJet, weight);
      histsBB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsBB.ptJet1->Fill(mva.ptJet1, weight);
      histsBB.ptJet2->Fill(mva.ptJet2, weight);
      histsBB.etaJet1->Fill(mva.etaJet1, weight);
      histsBB.etaJet2->Fill(mva.etaJet2, weight);
      histsBB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsBB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsBB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsBB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsBB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsBB.nJets->Fill(mva.nJets, weight);
      histsBB.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsBB.mDiMu->Fill(mva.mDiMu, weight);
      histsBB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsBB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsBB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsBB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsBB.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsBB.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsBB.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsBB.BDTHistVBF->Fill(bdtValVBF, weight);
        histsBB.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsBB.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (isBO)
    {
      histsBO.yDiMu->Fill(mva.yDiMu, weight);
      histsBO.ptDiMu->Fill(mva.ptDiMu, weight);
      histsBO.ptMu1->Fill(mva.ptMu1, weight);
      histsBO.ptMu2->Fill(mva.ptMu2, weight);
      histsBO.etaMu1->Fill(mva.etaMu1, weight);
      histsBO.etaMu2->Fill(mva.etaMu2, weight);
      histsBO.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsBO.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsBO.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsBO.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsBO.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsBO.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsBO.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsBO.nPU->Fill(nPU, weight);
      histsBO.nVtx->Fill(mva.nVtx, weight);
      histsBO.met->Fill(met.pt, weight);

      histsBO.mDiJet->Fill(mva.mDiJet, weight);
      histsBO.ptDiJet->Fill(mva.ptDiJet, weight);
      histsBO.ptJet1->Fill(mva.ptJet1, weight);
      histsBO.ptJet2->Fill(mva.ptJet2, weight);
      histsBO.etaJet1->Fill(mva.etaJet1, weight);
      histsBO.etaJet2->Fill(mva.etaJet2, weight);
      histsBO.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsBO.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsBO.deltaRJets->Fill(mva.deltaRJets, weight);
      histsBO.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsBO.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsBO.nJets->Fill(mva.nJets, weight);
      histsBO.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsBO.mDiMu->Fill(mva.mDiMu, weight);
      histsBO.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsBO.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsBO.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsBO.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsBO.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsBO.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsBO.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsBO.BDTHistVBF->Fill(bdtValVBF, weight);
        histsBO.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsBO.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (isBE)
    {
      histsBE.yDiMu->Fill(mva.yDiMu, weight);
      histsBE.ptDiMu->Fill(mva.ptDiMu, weight);
      histsBE.ptMu1->Fill(mva.ptMu1, weight);
      histsBE.ptMu2->Fill(mva.ptMu2, weight);
      histsBE.etaMu1->Fill(mva.etaMu1, weight);
      histsBE.etaMu2->Fill(mva.etaMu2, weight);
      histsBE.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsBE.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsBE.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsBE.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsBE.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsBE.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsBE.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsBE.nPU->Fill(nPU, weight);
      histsBE.nVtx->Fill(mva.nVtx, weight);
      histsBE.met->Fill(met.pt, weight);

      histsBE.mDiJet->Fill(mva.mDiJet, weight);
      histsBE.ptDiJet->Fill(mva.ptDiJet, weight);
      histsBE.ptJet1->Fill(mva.ptJet1, weight);
      histsBE.ptJet2->Fill(mva.ptJet2, weight);
      histsBE.etaJet1->Fill(mva.etaJet1, weight);
      histsBE.etaJet2->Fill(mva.etaJet2, weight);
      histsBE.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsBE.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsBE.deltaRJets->Fill(mva.deltaRJets, weight);
      histsBE.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsBE.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsBE.nJets->Fill(mva.nJets, weight);
      histsBE.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsBE.mDiMu->Fill(mva.mDiMu, weight);
      histsBE.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsBE.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsBE.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsBE.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsBE.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsBE.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsBE.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsBE.BDTHistVBF->Fill(bdtValVBF, weight);
        histsBE.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsBE.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (isOO)
    {
      histsOO.yDiMu->Fill(mva.yDiMu, weight);
      histsOO.ptDiMu->Fill(mva.ptDiMu, weight);
      histsOO.ptMu1->Fill(mva.ptMu1, weight);
      histsOO.ptMu2->Fill(mva.ptMu2, weight);
      histsOO.etaMu1->Fill(mva.etaMu1, weight);
      histsOO.etaMu2->Fill(mva.etaMu2, weight);
      histsOO.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsOO.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsOO.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsOO.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsOO.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsOO.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsOO.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsOO.nPU->Fill(nPU, weight);
      histsOO.nVtx->Fill(mva.nVtx, weight);
      histsOO.met->Fill(met.pt, weight);

      histsOO.mDiJet->Fill(mva.mDiJet, weight);
      histsOO.ptDiJet->Fill(mva.ptDiJet, weight);
      histsOO.ptJet1->Fill(mva.ptJet1, weight);
      histsOO.ptJet2->Fill(mva.ptJet2, weight);
      histsOO.etaJet1->Fill(mva.etaJet1, weight);
      histsOO.etaJet2->Fill(mva.etaJet2, weight);
      histsOO.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsOO.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsOO.deltaRJets->Fill(mva.deltaRJets, weight);
      histsOO.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsOO.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsOO.nJets->Fill(mva.nJets, weight);
      histsOO.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsOO.mDiMu->Fill(mva.mDiMu, weight);
      histsOO.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsOO.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsOO.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsOO.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsOO.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsOO.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsOO.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsOO.BDTHistVBF->Fill(bdtValVBF, weight);
        histsOO.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsOO.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (isOE)
    {
      histsOE.yDiMu->Fill(mva.yDiMu, weight);
      histsOE.ptDiMu->Fill(mva.ptDiMu, weight);
      histsOE.ptMu1->Fill(mva.ptMu1, weight);
      histsOE.ptMu2->Fill(mva.ptMu2, weight);
      histsOE.etaMu1->Fill(mva.etaMu1, weight);
      histsOE.etaMu2->Fill(mva.etaMu2, weight);
      histsOE.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsOE.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsOE.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsOE.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsOE.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsOE.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsOE.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsOE.nPU->Fill(nPU, weight);
      histsOE.nVtx->Fill(mva.nVtx, weight);
      histsOE.met->Fill(met.pt, weight);

      histsOE.mDiJet->Fill(mva.mDiJet, weight);
      histsOE.ptDiJet->Fill(mva.ptDiJet, weight);
      histsOE.ptJet1->Fill(mva.ptJet1, weight);
      histsOE.ptJet2->Fill(mva.ptJet2, weight);
      histsOE.etaJet1->Fill(mva.etaJet1, weight);
      histsOE.etaJet2->Fill(mva.etaJet2, weight);
      histsOE.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsOE.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsOE.deltaRJets->Fill(mva.deltaRJets, weight);
      histsOE.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsOE.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsOE.nJets->Fill(mva.nJets, weight);
      histsOE.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsOE.mDiMu->Fill(mva.mDiMu, weight);
      histsOE.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsOE.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsOE.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsOE.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsOE.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsOE.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsOE.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsOE.BDTHistVBF->Fill(bdtValVBF, weight);
        histsOE.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsOE.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (isEE)
    {
      histsEE.yDiMu->Fill(mva.yDiMu, weight);
      histsEE.ptDiMu->Fill(mva.ptDiMu, weight);
      histsEE.ptMu1->Fill(mva.ptMu1, weight);
      histsEE.ptMu2->Fill(mva.ptMu2, weight);
      histsEE.etaMu1->Fill(mva.etaMu1, weight);
      histsEE.etaMu2->Fill(mva.etaMu2, weight);
      histsEE.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsEE.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsEE.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsEE.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsEE.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsEE.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsEE.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsEE.nPU->Fill(nPU, weight);
      histsEE.nVtx->Fill(mva.nVtx, weight);
      histsEE.met->Fill(met.pt, weight);

      histsEE.mDiJet->Fill(mva.mDiJet, weight);
      histsEE.ptDiJet->Fill(mva.ptDiJet, weight);
      histsEE.ptJet1->Fill(mva.ptJet1, weight);
      histsEE.ptJet2->Fill(mva.ptJet2, weight);
      histsEE.etaJet1->Fill(mva.etaJet1, weight);
      histsEE.etaJet2->Fill(mva.etaJet2, weight);
      histsEE.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsEE.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsEE.deltaRJets->Fill(mva.deltaRJets, weight);
      histsEE.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsEE.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsEE.nJets->Fill(mva.nJets, weight);
      histsEE.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsEE.mDiMu->Fill(mva.mDiMu, weight);
      histsEE.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsEE.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsEE.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsEE.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsEE.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsEE.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsEE.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsEE.BDTHistVBF->Fill(bdtValVBF, weight);
        histsEE.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsEE.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (isNotBB)
    {
      histsNotBB.yDiMu->Fill(mva.yDiMu, weight);
      histsNotBB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsNotBB.ptMu1->Fill(mva.ptMu1, weight);
      histsNotBB.ptMu2->Fill(mva.ptMu2, weight);
      histsNotBB.etaMu1->Fill(mva.etaMu1, weight);
      histsNotBB.etaMu2->Fill(mva.etaMu2, weight);
      histsNotBB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsNotBB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsNotBB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsNotBB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsNotBB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsNotBB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsNotBB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsNotBB.nPU->Fill(nPU, weight);
      histsNotBB.nVtx->Fill(mva.nVtx, weight);
      histsNotBB.met->Fill(met.pt, weight);

      histsNotBB.mDiJet->Fill(mva.mDiJet, weight);
      histsNotBB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsNotBB.ptJet1->Fill(mva.ptJet1, weight);
      histsNotBB.ptJet2->Fill(mva.ptJet2, weight);
      histsNotBB.etaJet1->Fill(mva.etaJet1, weight);
      histsNotBB.etaJet2->Fill(mva.etaJet2, weight);
      histsNotBB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsNotBB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsNotBB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsNotBB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsNotBB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsNotBB.nJets->Fill(mva.nJets, weight);
      histsNotBB.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsNotBB.mDiMu->Fill(mva.mDiMu, weight);
      histsNotBB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsNotBB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsNotBB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsNotBB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsNotBB.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsNotBB.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsNotBB.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsNotBB.BDTHistVBF->Fill(bdtValVBF, weight);
        histsNotBB.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsNotBB.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
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
      histsVBFPresel.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsVBFPresel.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsVBFPresel.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsVBFPresel.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

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

    if (vbfPreselection && isBB)
    {
      histsVBFPreselBB.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFPreselBB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFPreselBB.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFPreselBB.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFPreselBB.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFPreselBB.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFPreselBB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFPreselBB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFPreselBB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFPreselBB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFPreselBB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFPreselBB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFPreselBB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFPreselBB.nPU->Fill(nPU, weight);
      histsVBFPreselBB.nVtx->Fill(mva.nVtx, weight);
      histsVBFPreselBB.met->Fill(met.pt, weight);

      histsVBFPreselBB.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFPreselBB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFPreselBB.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFPreselBB.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFPreselBB.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFPreselBB.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFPreselBB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFPreselBB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFPreselBB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFPreselBB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFPreselBB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFPreselBB.nJets->Fill(mva.nJets, weight);
      histsVBFPreselBB.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFPreselBB.mDiMu->Fill(mva.mDiMu, weight);
      histsVBFPreselBB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsVBFPreselBB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsVBFPreselBB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsVBFPreselBB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsVBFPreselBB.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsVBFPreselBB.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsVBFPreselBB.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsVBFPreselBB.BDTHistVBF->Fill(bdtValVBF, weight);
        histsVBFPreselBB.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsVBFPreselBB.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (vbfPreselection && isNotBB)
    {
      histsVBFPreselNotBB.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFPreselNotBB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFPreselNotBB.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFPreselNotBB.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFPreselNotBB.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFPreselNotBB.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFPreselNotBB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFPreselNotBB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFPreselNotBB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFPreselNotBB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFPreselNotBB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFPreselNotBB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFPreselNotBB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFPreselNotBB.nPU->Fill(nPU, weight);
      histsVBFPreselNotBB.nVtx->Fill(mva.nVtx, weight);
      histsVBFPreselNotBB.met->Fill(met.pt, weight);

      histsVBFPreselNotBB.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFPreselNotBB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFPreselNotBB.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFPreselNotBB.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFPreselNotBB.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFPreselNotBB.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFPreselNotBB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFPreselNotBB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFPreselNotBB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFPreselNotBB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFPreselNotBB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFPreselNotBB.nJets->Fill(mva.nJets, weight);
      histsVBFPreselNotBB.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFPreselNotBB.mDiMu->Fill(mva.mDiMu, weight);
      histsVBFPreselNotBB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsVBFPreselNotBB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsVBFPreselNotBB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsVBFPreselNotBB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsVBFPreselNotBB.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsVBFPreselNotBB.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsVBFPreselNotBB.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsVBFPreselNotBB.BDTHistVBF->Fill(bdtValVBF, weight);
        histsVBFPreselNotBB.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsVBFPreselNotBB.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
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
      histsIncPresel.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncPresel.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncPresel.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncPresel.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

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

    if (!vbfPreselection && isBB)
    {
      histsIncPreselBB.yDiMu->Fill(mva.yDiMu, weight);
      histsIncPreselBB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncPreselBB.ptMu1->Fill(mva.ptMu1, weight);
      histsIncPreselBB.ptMu2->Fill(mva.ptMu2, weight);
      histsIncPreselBB.etaMu1->Fill(mva.etaMu1, weight);
      histsIncPreselBB.etaMu2->Fill(mva.etaMu2, weight);
      histsIncPreselBB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncPreselBB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncPreselBB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncPreselBB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncPreselBB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncPreselBB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncPreselBB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncPreselBB.nPU->Fill(nPU, weight);
      histsIncPreselBB.nVtx->Fill(mva.nVtx, weight);
      histsIncPreselBB.met->Fill(met.pt, weight);

      histsIncPreselBB.mDiJet->Fill(mva.mDiJet, weight);
      histsIncPreselBB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncPreselBB.ptJet1->Fill(mva.ptJet1, weight);
      histsIncPreselBB.ptJet2->Fill(mva.ptJet2, weight);
      histsIncPreselBB.etaJet1->Fill(mva.etaJet1, weight);
      histsIncPreselBB.etaJet2->Fill(mva.etaJet2, weight);
      histsIncPreselBB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncPreselBB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncPreselBB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncPreselBB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncPreselBB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncPreselBB.nJets->Fill(mva.nJets, weight);
      histsIncPreselBB.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncPreselBB.mDiMu->Fill(mva.mDiMu, weight);
      histsIncPreselBB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncPreselBB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncPreselBB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncPreselBB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsIncPreselBB.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsIncPreselBB.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsIncPreselBB.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsIncPreselBB.BDTHistVBF->Fill(bdtValVBF, weight);
        histsIncPreselBB.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsIncPreselBB.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (!vbfPreselection && isBO)
    {
      histsIncPreselBO.yDiMu->Fill(mva.yDiMu, weight);
      histsIncPreselBO.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncPreselBO.ptMu1->Fill(mva.ptMu1, weight);
      histsIncPreselBO.ptMu2->Fill(mva.ptMu2, weight);
      histsIncPreselBO.etaMu1->Fill(mva.etaMu1, weight);
      histsIncPreselBO.etaMu2->Fill(mva.etaMu2, weight);
      histsIncPreselBO.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncPreselBO.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncPreselBO.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncPreselBO.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncPreselBO.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncPreselBO.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncPreselBO.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncPreselBO.nPU->Fill(nPU, weight);
      histsIncPreselBO.nVtx->Fill(mva.nVtx, weight);
      histsIncPreselBO.met->Fill(met.pt, weight);

      histsIncPreselBO.mDiJet->Fill(mva.mDiJet, weight);
      histsIncPreselBO.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncPreselBO.ptJet1->Fill(mva.ptJet1, weight);
      histsIncPreselBO.ptJet2->Fill(mva.ptJet2, weight);
      histsIncPreselBO.etaJet1->Fill(mva.etaJet1, weight);
      histsIncPreselBO.etaJet2->Fill(mva.etaJet2, weight);
      histsIncPreselBO.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncPreselBO.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncPreselBO.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncPreselBO.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncPreselBO.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncPreselBO.nJets->Fill(mva.nJets, weight);
      histsIncPreselBO.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncPreselBO.mDiMu->Fill(mva.mDiMu, weight);
      histsIncPreselBO.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncPreselBO.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncPreselBO.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncPreselBO.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsIncPreselBO.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsIncPreselBO.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsIncPreselBO.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsIncPreselBO.BDTHistVBF->Fill(bdtValVBF, weight);
        histsIncPreselBO.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsIncPreselBO.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (!vbfPreselection && isBE)
    {
      histsIncPreselBE.yDiMu->Fill(mva.yDiMu, weight);
      histsIncPreselBE.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncPreselBE.ptMu1->Fill(mva.ptMu1, weight);
      histsIncPreselBE.ptMu2->Fill(mva.ptMu2, weight);
      histsIncPreselBE.etaMu1->Fill(mva.etaMu1, weight);
      histsIncPreselBE.etaMu2->Fill(mva.etaMu2, weight);
      histsIncPreselBE.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncPreselBE.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncPreselBE.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncPreselBE.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncPreselBE.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncPreselBE.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncPreselBE.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncPreselBE.nPU->Fill(nPU, weight);
      histsIncPreselBE.nVtx->Fill(mva.nVtx, weight);
      histsIncPreselBE.met->Fill(met.pt, weight);

      histsIncPreselBE.mDiJet->Fill(mva.mDiJet, weight);
      histsIncPreselBE.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncPreselBE.ptJet1->Fill(mva.ptJet1, weight);
      histsIncPreselBE.ptJet2->Fill(mva.ptJet2, weight);
      histsIncPreselBE.etaJet1->Fill(mva.etaJet1, weight);
      histsIncPreselBE.etaJet2->Fill(mva.etaJet2, weight);
      histsIncPreselBE.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncPreselBE.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncPreselBE.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncPreselBE.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncPreselBE.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncPreselBE.nJets->Fill(mva.nJets, weight);
      histsIncPreselBE.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncPreselBE.mDiMu->Fill(mva.mDiMu, weight);
      histsIncPreselBE.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncPreselBE.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncPreselBE.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncPreselBE.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsIncPreselBE.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsIncPreselBE.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsIncPreselBE.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsIncPreselBE.BDTHistVBF->Fill(bdtValVBF, weight);
        histsIncPreselBE.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsIncPreselBE.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (!vbfPreselection && isOO)
    {
      histsIncPreselOO.yDiMu->Fill(mva.yDiMu, weight);
      histsIncPreselOO.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncPreselOO.ptMu1->Fill(mva.ptMu1, weight);
      histsIncPreselOO.ptMu2->Fill(mva.ptMu2, weight);
      histsIncPreselOO.etaMu1->Fill(mva.etaMu1, weight);
      histsIncPreselOO.etaMu2->Fill(mva.etaMu2, weight);
      histsIncPreselOO.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncPreselOO.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncPreselOO.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncPreselOO.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncPreselOO.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncPreselOO.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncPreselOO.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncPreselOO.nPU->Fill(nPU, weight);
      histsIncPreselOO.nVtx->Fill(mva.nVtx, weight);
      histsIncPreselOO.met->Fill(met.pt, weight);

      histsIncPreselOO.mDiJet->Fill(mva.mDiJet, weight);
      histsIncPreselOO.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncPreselOO.ptJet1->Fill(mva.ptJet1, weight);
      histsIncPreselOO.ptJet2->Fill(mva.ptJet2, weight);
      histsIncPreselOO.etaJet1->Fill(mva.etaJet1, weight);
      histsIncPreselOO.etaJet2->Fill(mva.etaJet2, weight);
      histsIncPreselOO.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncPreselOO.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncPreselOO.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncPreselOO.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncPreselOO.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncPreselOO.nJets->Fill(mva.nJets, weight);
      histsIncPreselOO.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncPreselOO.mDiMu->Fill(mva.mDiMu, weight);
      histsIncPreselOO.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncPreselOO.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncPreselOO.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncPreselOO.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsIncPreselOO.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsIncPreselOO.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsIncPreselOO.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsIncPreselOO.BDTHistVBF->Fill(bdtValVBF, weight);
        histsIncPreselOO.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsIncPreselOO.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (!vbfPreselection && isOE)
    {
      histsIncPreselOE.yDiMu->Fill(mva.yDiMu, weight);
      histsIncPreselOE.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncPreselOE.ptMu1->Fill(mva.ptMu1, weight);
      histsIncPreselOE.ptMu2->Fill(mva.ptMu2, weight);
      histsIncPreselOE.etaMu1->Fill(mva.etaMu1, weight);
      histsIncPreselOE.etaMu2->Fill(mva.etaMu2, weight);
      histsIncPreselOE.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncPreselOE.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncPreselOE.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncPreselOE.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncPreselOE.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncPreselOE.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncPreselOE.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncPreselOE.nPU->Fill(nPU, weight);
      histsIncPreselOE.nVtx->Fill(mva.nVtx, weight);
      histsIncPreselOE.met->Fill(met.pt, weight);

      histsIncPreselOE.mDiJet->Fill(mva.mDiJet, weight);
      histsIncPreselOE.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncPreselOE.ptJet1->Fill(mva.ptJet1, weight);
      histsIncPreselOE.ptJet2->Fill(mva.ptJet2, weight);
      histsIncPreselOE.etaJet1->Fill(mva.etaJet1, weight);
      histsIncPreselOE.etaJet2->Fill(mva.etaJet2, weight);
      histsIncPreselOE.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncPreselOE.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncPreselOE.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncPreselOE.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncPreselOE.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncPreselOE.nJets->Fill(mva.nJets, weight);
      histsIncPreselOE.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncPreselOE.mDiMu->Fill(mva.mDiMu, weight);
      histsIncPreselOE.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncPreselOE.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncPreselOE.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncPreselOE.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsIncPreselOE.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsIncPreselOE.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsIncPreselOE.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsIncPreselOE.BDTHistVBF->Fill(bdtValVBF, weight);
        histsIncPreselOE.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsIncPreselOE.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (!vbfPreselection && isEE)
    {
      histsIncPreselEE.yDiMu->Fill(mva.yDiMu, weight);
      histsIncPreselEE.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncPreselEE.ptMu1->Fill(mva.ptMu1, weight);
      histsIncPreselEE.ptMu2->Fill(mva.ptMu2, weight);
      histsIncPreselEE.etaMu1->Fill(mva.etaMu1, weight);
      histsIncPreselEE.etaMu2->Fill(mva.etaMu2, weight);
      histsIncPreselEE.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncPreselEE.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncPreselEE.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncPreselEE.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncPreselEE.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncPreselEE.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncPreselEE.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncPreselEE.nPU->Fill(nPU, weight);
      histsIncPreselEE.nVtx->Fill(mva.nVtx, weight);
      histsIncPreselEE.met->Fill(met.pt, weight);

      histsIncPreselEE.mDiJet->Fill(mva.mDiJet, weight);
      histsIncPreselEE.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncPreselEE.ptJet1->Fill(mva.ptJet1, weight);
      histsIncPreselEE.ptJet2->Fill(mva.ptJet2, weight);
      histsIncPreselEE.etaJet1->Fill(mva.etaJet1, weight);
      histsIncPreselEE.etaJet2->Fill(mva.etaJet2, weight);
      histsIncPreselEE.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncPreselEE.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncPreselEE.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncPreselEE.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncPreselEE.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncPreselEE.nJets->Fill(mva.nJets, weight);
      histsIncPreselEE.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncPreselEE.mDiMu->Fill(mva.mDiMu, weight);
      histsIncPreselEE.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncPreselEE.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncPreselEE.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncPreselEE.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsIncPreselEE.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsIncPreselEE.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsIncPreselEE.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsIncPreselEE.BDTHistVBF->Fill(bdtValVBF, weight);
        histsIncPreselEE.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsIncPreselEE.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (!vbfPreselection && isNotBB)
    {
      histsIncPreselNotBB.yDiMu->Fill(mva.yDiMu, weight);
      histsIncPreselNotBB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncPreselNotBB.ptMu1->Fill(mva.ptMu1, weight);
      histsIncPreselNotBB.ptMu2->Fill(mva.ptMu2, weight);
      histsIncPreselNotBB.etaMu1->Fill(mva.etaMu1, weight);
      histsIncPreselNotBB.etaMu2->Fill(mva.etaMu2, weight);
      histsIncPreselNotBB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncPreselNotBB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncPreselNotBB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncPreselNotBB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncPreselNotBB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncPreselNotBB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncPreselNotBB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncPreselNotBB.nPU->Fill(nPU, weight);
      histsIncPreselNotBB.nVtx->Fill(mva.nVtx, weight);
      histsIncPreselNotBB.met->Fill(met.pt, weight);

      histsIncPreselNotBB.mDiJet->Fill(mva.mDiJet, weight);
      histsIncPreselNotBB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncPreselNotBB.ptJet1->Fill(mva.ptJet1, weight);
      histsIncPreselNotBB.ptJet2->Fill(mva.ptJet2, weight);
      histsIncPreselNotBB.etaJet1->Fill(mva.etaJet1, weight);
      histsIncPreselNotBB.etaJet2->Fill(mva.etaJet2, weight);
      histsIncPreselNotBB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncPreselNotBB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncPreselNotBB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncPreselNotBB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncPreselNotBB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncPreselNotBB.nJets->Fill(mva.nJets, weight);
      histsIncPreselNotBB.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncPreselNotBB.mDiMu->Fill(mva.mDiMu, weight);
      histsIncPreselNotBB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncPreselNotBB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncPreselNotBB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncPreselNotBB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsIncPreselNotBB.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsIncPreselNotBB.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsIncPreselNotBB.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsIncPreselNotBB.BDTHistVBF->Fill(bdtValVBF, weight);
        histsIncPreselNotBB.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsIncPreselNotBB.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (!vbfPreselection && mva.getMVAPassBDTCut(cfgNameInc))
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncBDTSig80.mDiMu->Fill(mva.mDiMu, weight);
      histsIncBDTSig80.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncBDTSig80.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncBDTSig80.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncBDTSig80.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

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

    if (!vbfPreselection && mva.getMVAPassBDTCut(cfgNameInc) && isBB)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncBDTSig80BB.mDiMu->Fill(mva.mDiMu, weight);
      histsIncBDTSig80BB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncBDTSig80BB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncBDTSig80BB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncBDTSig80BB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

#ifdef BLIND
      }
#endif
      histsIncBDTSig80BB.yDiMu->Fill(mva.yDiMu, weight);
      histsIncBDTSig80BB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncBDTSig80BB.ptMu1->Fill(mva.ptMu1, weight);
      histsIncBDTSig80BB.ptMu2->Fill(mva.ptMu2, weight);
      histsIncBDTSig80BB.etaMu1->Fill(mva.etaMu1, weight);
      histsIncBDTSig80BB.etaMu2->Fill(mva.etaMu2, weight);
      histsIncBDTSig80BB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncBDTSig80BB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncBDTSig80BB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncBDTSig80BB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncBDTSig80BB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncBDTSig80BB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncBDTSig80BB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncBDTSig80BB.nPU->Fill(nPU, weight);
      histsIncBDTSig80BB.nVtx->Fill(mva.nVtx, weight);
      histsIncBDTSig80BB.met->Fill(met.pt, weight);

      histsIncBDTSig80BB.mDiJet->Fill(mva.mDiJet, weight);
      histsIncBDTSig80BB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncBDTSig80BB.ptJet1->Fill(mva.ptJet1, weight);
      histsIncBDTSig80BB.ptJet2->Fill(mva.ptJet2, weight);
      histsIncBDTSig80BB.etaJet1->Fill(mva.etaJet1, weight);
      histsIncBDTSig80BB.etaJet2->Fill(mva.etaJet2, weight);
      histsIncBDTSig80BB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncBDTSig80BB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncBDTSig80BB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncBDTSig80BB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncBDTSig80BB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncBDTSig80BB.nJets->Fill(mva.nJets, weight);
      histsIncBDTSig80BB.ht->Fill(mva.ht, weight);
      histsIncBDTSig80BB.BDTHistMuonOnly->Fill(bdtValInc, weight);
    }

    if (!vbfPreselection && mva.getMVAPassBDTCut(cfgNameInc) && isBO)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncBDTSig80BO.mDiMu->Fill(mva.mDiMu, weight);
      histsIncBDTSig80BO.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncBDTSig80BO.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncBDTSig80BO.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncBDTSig80BO.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

#ifdef BLIND
      }
#endif
      histsIncBDTSig80BO.yDiMu->Fill(mva.yDiMu, weight);
      histsIncBDTSig80BO.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncBDTSig80BO.ptMu1->Fill(mva.ptMu1, weight);
      histsIncBDTSig80BO.ptMu2->Fill(mva.ptMu2, weight);
      histsIncBDTSig80BO.etaMu1->Fill(mva.etaMu1, weight);
      histsIncBDTSig80BO.etaMu2->Fill(mva.etaMu2, weight);
      histsIncBDTSig80BO.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncBDTSig80BO.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncBDTSig80BO.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncBDTSig80BO.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncBDTSig80BO.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncBDTSig80BO.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncBDTSig80BO.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncBDTSig80BO.nPU->Fill(nPU, weight);
      histsIncBDTSig80BO.nVtx->Fill(mva.nVtx, weight);
      histsIncBDTSig80BO.met->Fill(met.pt, weight);

      histsIncBDTSig80BO.mDiJet->Fill(mva.mDiJet, weight);
      histsIncBDTSig80BO.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncBDTSig80BO.ptJet1->Fill(mva.ptJet1, weight);
      histsIncBDTSig80BO.ptJet2->Fill(mva.ptJet2, weight);
      histsIncBDTSig80BO.etaJet1->Fill(mva.etaJet1, weight);
      histsIncBDTSig80BO.etaJet2->Fill(mva.etaJet2, weight);
      histsIncBDTSig80BO.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncBDTSig80BO.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncBDTSig80BO.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncBDTSig80BO.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncBDTSig80BO.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncBDTSig80BO.nJets->Fill(mva.nJets, weight);
      histsIncBDTSig80BO.ht->Fill(mva.ht, weight);
      histsIncBDTSig80BO.BDTHistMuonOnly->Fill(bdtValInc, weight);
    }

    if (!vbfPreselection && mva.getMVAPassBDTCut(cfgNameInc) && isBE)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncBDTSig80BE.mDiMu->Fill(mva.mDiMu, weight);
      histsIncBDTSig80BE.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncBDTSig80BE.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncBDTSig80BE.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncBDTSig80BE.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

#ifdef BLIND
      }
#endif
      histsIncBDTSig80BE.yDiMu->Fill(mva.yDiMu, weight);
      histsIncBDTSig80BE.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncBDTSig80BE.ptMu1->Fill(mva.ptMu1, weight);
      histsIncBDTSig80BE.ptMu2->Fill(mva.ptMu2, weight);
      histsIncBDTSig80BE.etaMu1->Fill(mva.etaMu1, weight);
      histsIncBDTSig80BE.etaMu2->Fill(mva.etaMu2, weight);
      histsIncBDTSig80BE.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncBDTSig80BE.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncBDTSig80BE.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncBDTSig80BE.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncBDTSig80BE.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncBDTSig80BE.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncBDTSig80BE.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncBDTSig80BE.nPU->Fill(nPU, weight);
      histsIncBDTSig80BE.nVtx->Fill(mva.nVtx, weight);
      histsIncBDTSig80BE.met->Fill(met.pt, weight);

      histsIncBDTSig80BE.mDiJet->Fill(mva.mDiJet, weight);
      histsIncBDTSig80BE.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncBDTSig80BE.ptJet1->Fill(mva.ptJet1, weight);
      histsIncBDTSig80BE.ptJet2->Fill(mva.ptJet2, weight);
      histsIncBDTSig80BE.etaJet1->Fill(mva.etaJet1, weight);
      histsIncBDTSig80BE.etaJet2->Fill(mva.etaJet2, weight);
      histsIncBDTSig80BE.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncBDTSig80BE.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncBDTSig80BE.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncBDTSig80BE.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncBDTSig80BE.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncBDTSig80BE.nJets->Fill(mva.nJets, weight);
      histsIncBDTSig80BE.ht->Fill(mva.ht, weight);
      histsIncBDTSig80BE.BDTHistMuonOnly->Fill(bdtValInc, weight);
    }

    if (!vbfPreselection && mva.getMVAPassBDTCut(cfgNameInc) && isOO)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncBDTSig80OO.mDiMu->Fill(mva.mDiMu, weight);
      histsIncBDTSig80OO.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncBDTSig80OO.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncBDTSig80OO.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncBDTSig80OO.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

#ifdef BLIND
      }
#endif
      histsIncBDTSig80OO.yDiMu->Fill(mva.yDiMu, weight);
      histsIncBDTSig80OO.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncBDTSig80OO.ptMu1->Fill(mva.ptMu1, weight);
      histsIncBDTSig80OO.ptMu2->Fill(mva.ptMu2, weight);
      histsIncBDTSig80OO.etaMu1->Fill(mva.etaMu1, weight);
      histsIncBDTSig80OO.etaMu2->Fill(mva.etaMu2, weight);
      histsIncBDTSig80OO.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncBDTSig80OO.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncBDTSig80OO.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncBDTSig80OO.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncBDTSig80OO.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncBDTSig80OO.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncBDTSig80OO.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncBDTSig80OO.nPU->Fill(nPU, weight);
      histsIncBDTSig80OO.nVtx->Fill(mva.nVtx, weight);
      histsIncBDTSig80OO.met->Fill(met.pt, weight);

      histsIncBDTSig80OO.mDiJet->Fill(mva.mDiJet, weight);
      histsIncBDTSig80OO.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncBDTSig80OO.ptJet1->Fill(mva.ptJet1, weight);
      histsIncBDTSig80OO.ptJet2->Fill(mva.ptJet2, weight);
      histsIncBDTSig80OO.etaJet1->Fill(mva.etaJet1, weight);
      histsIncBDTSig80OO.etaJet2->Fill(mva.etaJet2, weight);
      histsIncBDTSig80OO.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncBDTSig80OO.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncBDTSig80OO.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncBDTSig80OO.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncBDTSig80OO.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncBDTSig80OO.nJets->Fill(mva.nJets, weight);
      histsIncBDTSig80OO.ht->Fill(mva.ht, weight);
      histsIncBDTSig80OO.BDTHistMuonOnly->Fill(bdtValInc, weight);
    }

    if (!vbfPreselection && mva.getMVAPassBDTCut(cfgNameInc) && isOE)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncBDTSig80OE.mDiMu->Fill(mva.mDiMu, weight);
      histsIncBDTSig80OE.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncBDTSig80OE.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncBDTSig80OE.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncBDTSig80OE.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

#ifdef BLIND
      }
#endif
      histsIncBDTSig80OE.yDiMu->Fill(mva.yDiMu, weight);
      histsIncBDTSig80OE.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncBDTSig80OE.ptMu1->Fill(mva.ptMu1, weight);
      histsIncBDTSig80OE.ptMu2->Fill(mva.ptMu2, weight);
      histsIncBDTSig80OE.etaMu1->Fill(mva.etaMu1, weight);
      histsIncBDTSig80OE.etaMu2->Fill(mva.etaMu2, weight);
      histsIncBDTSig80OE.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncBDTSig80OE.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncBDTSig80OE.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncBDTSig80OE.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncBDTSig80OE.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncBDTSig80OE.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncBDTSig80OE.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncBDTSig80OE.nPU->Fill(nPU, weight);
      histsIncBDTSig80OE.nVtx->Fill(mva.nVtx, weight);
      histsIncBDTSig80OE.met->Fill(met.pt, weight);

      histsIncBDTSig80OE.mDiJet->Fill(mva.mDiJet, weight);
      histsIncBDTSig80OE.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncBDTSig80OE.ptJet1->Fill(mva.ptJet1, weight);
      histsIncBDTSig80OE.ptJet2->Fill(mva.ptJet2, weight);
      histsIncBDTSig80OE.etaJet1->Fill(mva.etaJet1, weight);
      histsIncBDTSig80OE.etaJet2->Fill(mva.etaJet2, weight);
      histsIncBDTSig80OE.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncBDTSig80OE.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncBDTSig80OE.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncBDTSig80OE.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncBDTSig80OE.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncBDTSig80OE.nJets->Fill(mva.nJets, weight);
      histsIncBDTSig80OE.ht->Fill(mva.ht, weight);
      histsIncBDTSig80OE.BDTHistMuonOnly->Fill(bdtValInc, weight);
    }

    if (!vbfPreselection && mva.getMVAPassBDTCut(cfgNameInc) && isEE)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncBDTSig80EE.mDiMu->Fill(mva.mDiMu, weight);
      histsIncBDTSig80EE.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncBDTSig80EE.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncBDTSig80EE.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncBDTSig80EE.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

#ifdef BLIND
      }
#endif
      histsIncBDTSig80EE.yDiMu->Fill(mva.yDiMu, weight);
      histsIncBDTSig80EE.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncBDTSig80EE.ptMu1->Fill(mva.ptMu1, weight);
      histsIncBDTSig80EE.ptMu2->Fill(mva.ptMu2, weight);
      histsIncBDTSig80EE.etaMu1->Fill(mva.etaMu1, weight);
      histsIncBDTSig80EE.etaMu2->Fill(mva.etaMu2, weight);
      histsIncBDTSig80EE.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncBDTSig80EE.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncBDTSig80EE.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncBDTSig80EE.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncBDTSig80EE.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncBDTSig80EE.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncBDTSig80EE.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncBDTSig80EE.nPU->Fill(nPU, weight);
      histsIncBDTSig80EE.nVtx->Fill(mva.nVtx, weight);
      histsIncBDTSig80EE.met->Fill(met.pt, weight);

      histsIncBDTSig80EE.mDiJet->Fill(mva.mDiJet, weight);
      histsIncBDTSig80EE.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncBDTSig80EE.ptJet1->Fill(mva.ptJet1, weight);
      histsIncBDTSig80EE.ptJet2->Fill(mva.ptJet2, weight);
      histsIncBDTSig80EE.etaJet1->Fill(mva.etaJet1, weight);
      histsIncBDTSig80EE.etaJet2->Fill(mva.etaJet2, weight);
      histsIncBDTSig80EE.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncBDTSig80EE.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncBDTSig80EE.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncBDTSig80EE.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncBDTSig80EE.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncBDTSig80EE.nJets->Fill(mva.nJets, weight);
      histsIncBDTSig80EE.ht->Fill(mva.ht, weight);
      histsIncBDTSig80EE.BDTHistMuonOnly->Fill(bdtValInc, weight);
    }

    if (!vbfPreselection && mva.getMVAPassBDTCut(cfgNameInc) && isNotBB)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncBDTSig80NotBB.mDiMu->Fill(mva.mDiMu, weight);
      histsIncBDTSig80NotBB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncBDTSig80NotBB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncBDTSig80NotBB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncBDTSig80NotBB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

#ifdef BLIND
      }
#endif
      histsIncBDTSig80NotBB.yDiMu->Fill(mva.yDiMu, weight);
      histsIncBDTSig80NotBB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncBDTSig80NotBB.ptMu1->Fill(mva.ptMu1, weight);
      histsIncBDTSig80NotBB.ptMu2->Fill(mva.ptMu2, weight);
      histsIncBDTSig80NotBB.etaMu1->Fill(mva.etaMu1, weight);
      histsIncBDTSig80NotBB.etaMu2->Fill(mva.etaMu2, weight);
      histsIncBDTSig80NotBB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncBDTSig80NotBB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncBDTSig80NotBB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncBDTSig80NotBB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncBDTSig80NotBB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncBDTSig80NotBB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncBDTSig80NotBB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncBDTSig80NotBB.nPU->Fill(nPU, weight);
      histsIncBDTSig80NotBB.nVtx->Fill(mva.nVtx, weight);
      histsIncBDTSig80NotBB.met->Fill(met.pt, weight);

      histsIncBDTSig80NotBB.mDiJet->Fill(mva.mDiJet, weight);
      histsIncBDTSig80NotBB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncBDTSig80NotBB.ptJet1->Fill(mva.ptJet1, weight);
      histsIncBDTSig80NotBB.ptJet2->Fill(mva.ptJet2, weight);
      histsIncBDTSig80NotBB.etaJet1->Fill(mva.etaJet1, weight);
      histsIncBDTSig80NotBB.etaJet2->Fill(mva.etaJet2, weight);
      histsIncBDTSig80NotBB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncBDTSig80NotBB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncBDTSig80NotBB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncBDTSig80NotBB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncBDTSig80NotBB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncBDTSig80NotBB.nJets->Fill(mva.nJets, weight);
      histsIncBDTSig80NotBB.ht->Fill(mva.ht, weight);
      histsIncBDTSig80NotBB.BDTHistMuonOnly->Fill(bdtValInc, weight);
    }

    if (vbfPreselection && mva.getMVAPassBDTCut(cfgNameVBF))
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFBDTSig80.mDiMu->Fill(mva.mDiMu, weight);
      histsVBFBDTSig80.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsVBFBDTSig80.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsVBFBDTSig80.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsVBFBDTSig80.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

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

    if (vbfPreselection && mva.getMVAPassBDTCut(cfgNameVBF) && isBB)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFBDTSig80BB.mDiMu->Fill(mva.mDiMu, weight);
      histsVBFBDTSig80BB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsVBFBDTSig80BB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsVBFBDTSig80BB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsVBFBDTSig80BB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

#ifdef BLIND
      }
#endif
      histsVBFBDTSig80BB.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFBDTSig80BB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFBDTSig80BB.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFBDTSig80BB.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFBDTSig80BB.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFBDTSig80BB.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFBDTSig80BB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFBDTSig80BB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFBDTSig80BB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFBDTSig80BB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFBDTSig80BB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFBDTSig80BB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFBDTSig80BB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFBDTSig80BB.nPU->Fill(nPU, weight);
      histsVBFBDTSig80BB.nVtx->Fill(mva.nVtx, weight);
      histsVBFBDTSig80BB.met->Fill(met.pt, weight);

      histsVBFBDTSig80BB.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFBDTSig80BB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFBDTSig80BB.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFBDTSig80BB.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFBDTSig80BB.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFBDTSig80BB.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFBDTSig80BB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFBDTSig80BB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFBDTSig80BB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFBDTSig80BB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFBDTSig80BB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFBDTSig80BB.nJets->Fill(mva.nJets, weight);
      histsVBFBDTSig80BB.ht->Fill(mva.ht, weight);
      histsVBFBDTSig80BB.BDTHistVBF->Fill(bdtValVBF, weight);
    }

    if (vbfPreselection && mva.getMVAPassBDTCut(cfgNameVBF) && isNotBB)
    {
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFBDTSig80NotBB.mDiMu->Fill(mva.mDiMu, weight);
      histsVBFBDTSig80NotBB.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsVBFBDTSig80NotBB.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsVBFBDTSig80NotBB.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsVBFBDTSig80NotBB.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

#ifdef BLIND
      }
#endif
      histsVBFBDTSig80NotBB.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFBDTSig80NotBB.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFBDTSig80NotBB.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFBDTSig80NotBB.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFBDTSig80NotBB.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFBDTSig80NotBB.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFBDTSig80NotBB.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFBDTSig80NotBB.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFBDTSig80NotBB.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFBDTSig80NotBB.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFBDTSig80NotBB.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFBDTSig80NotBB.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFBDTSig80NotBB.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFBDTSig80NotBB.nPU->Fill(nPU, weight);
      histsVBFBDTSig80NotBB.nVtx->Fill(mva.nVtx, weight);
      histsVBFBDTSig80NotBB.met->Fill(met.pt, weight);

      histsVBFBDTSig80NotBB.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFBDTSig80NotBB.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFBDTSig80NotBB.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFBDTSig80NotBB.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFBDTSig80NotBB.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFBDTSig80NotBB.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFBDTSig80NotBB.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFBDTSig80NotBB.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFBDTSig80NotBB.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFBDTSig80NotBB.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFBDTSig80NotBB.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFBDTSig80NotBB.nJets->Fill(mva.nJets, weight);
      histsVBFBDTSig80NotBB.ht->Fill(mva.ht, weight);
      histsVBFBDTSig80NotBB.BDTHistVBF->Fill(bdtValVBF, weight);
    }

    if (!vbfPreselection && mva.ptDiMu < 20.0)
    {
      histsIncPreselDiMuPtL20.yDiMu->Fill(mva.yDiMu, weight);
      histsIncPreselDiMuPtL20.ptDiMu->Fill(mva.ptDiMu, weight);
      histsIncPreselDiMuPtL20.ptMu1->Fill(mva.ptMu1, weight);
      histsIncPreselDiMuPtL20.ptMu2->Fill(mva.ptMu2, weight);
      histsIncPreselDiMuPtL20.etaMu1->Fill(mva.etaMu1, weight);
      histsIncPreselDiMuPtL20.etaMu2->Fill(mva.etaMu2, weight);
      histsIncPreselDiMuPtL20.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsIncPreselDiMuPtL20.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsIncPreselDiMuPtL20.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsIncPreselDiMuPtL20.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsIncPreselDiMuPtL20.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsIncPreselDiMuPtL20.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsIncPreselDiMuPtL20.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsIncPreselDiMuPtL20.nPU->Fill(nPU, weight);
      histsIncPreselDiMuPtL20.nVtx->Fill(mva.nVtx, weight);
      histsIncPreselDiMuPtL20.met->Fill(met.pt, weight);

      histsIncPreselDiMuPtL20.mDiJet->Fill(mva.mDiJet, weight);
      histsIncPreselDiMuPtL20.ptDiJet->Fill(mva.ptDiJet, weight);
      histsIncPreselDiMuPtL20.ptJet1->Fill(mva.ptJet1, weight);
      histsIncPreselDiMuPtL20.ptJet2->Fill(mva.ptJet2, weight);
      histsIncPreselDiMuPtL20.etaJet1->Fill(mva.etaJet1, weight);
      histsIncPreselDiMuPtL20.etaJet2->Fill(mva.etaJet2, weight);
      histsIncPreselDiMuPtL20.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsIncPreselDiMuPtL20.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsIncPreselDiMuPtL20.deltaRJets->Fill(mva.deltaRJets, weight);
      histsIncPreselDiMuPtL20.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsIncPreselDiMuPtL20.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsIncPreselDiMuPtL20.nJets->Fill(mva.nJets, weight);
      histsIncPreselDiMuPtL20.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsIncPreselDiMuPtL20.mDiMu->Fill(mva.mDiMu, weight);
      histsIncPreselDiMuPtL20.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsIncPreselDiMuPtL20.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsIncPreselDiMuPtL20.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsIncPreselDiMuPtL20.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsIncPreselDiMuPtL20.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsIncPreselDiMuPtL20.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsIncPreselDiMuPtL20.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsIncPreselDiMuPtL20.BDTHistVBF->Fill(bdtValVBF, weight);
        histsIncPreselDiMuPtL20.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsIncPreselDiMuPtL20.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    if (vbfPreselection && mva.ptDiMu < 20.0)
    {
      histsVBFPreselDiMuPtL20.yDiMu->Fill(mva.yDiMu, weight);
      histsVBFPreselDiMuPtL20.ptDiMu->Fill(mva.ptDiMu, weight);
      histsVBFPreselDiMuPtL20.ptMu1->Fill(mva.ptMu1, weight);
      histsVBFPreselDiMuPtL20.ptMu2->Fill(mva.ptMu2, weight);
      histsVBFPreselDiMuPtL20.etaMu1->Fill(mva.etaMu1, weight);
      histsVBFPreselDiMuPtL20.etaMu2->Fill(mva.etaMu2, weight);
      histsVBFPreselDiMuPtL20.cosThetaStar->Fill(mva.cosThetaStar, weight);
      histsVBFPreselDiMuPtL20.cosThetaStarCS->Fill(mva.cosThetaStarCS, weight);
      histsVBFPreselDiMuPtL20.deltaPhiMuons->Fill(mva.deltaPhiMuons, weight);
      histsVBFPreselDiMuPtL20.deltaEtaMuons->Fill(mva.deltaEtaMuons, weight);
      histsVBFPreselDiMuPtL20.deltaRMuons->Fill(mva.deltaRMuons, weight);
      histsVBFPreselDiMuPtL20.relIsoMu1->Fill(mva.relIsoMu1, weight);
      histsVBFPreselDiMuPtL20.relIsoMu2->Fill(mva.relIsoMu2, weight);
      histsVBFPreselDiMuPtL20.nPU->Fill(nPU, weight);
      histsVBFPreselDiMuPtL20.nVtx->Fill(mva.nVtx, weight);
      histsVBFPreselDiMuPtL20.met->Fill(met.pt, weight);

      histsVBFPreselDiMuPtL20.mDiJet->Fill(mva.mDiJet, weight);
      histsVBFPreselDiMuPtL20.ptDiJet->Fill(mva.ptDiJet, weight);
      histsVBFPreselDiMuPtL20.ptJet1->Fill(mva.ptJet1, weight);
      histsVBFPreselDiMuPtL20.ptJet2->Fill(mva.ptJet2, weight);
      histsVBFPreselDiMuPtL20.etaJet1->Fill(mva.etaJet1, weight);
      histsVBFPreselDiMuPtL20.etaJet2->Fill(mva.etaJet2, weight);
      histsVBFPreselDiMuPtL20.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      histsVBFPreselDiMuPtL20.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      histsVBFPreselDiMuPtL20.deltaRJets->Fill(mva.deltaRJets, weight);
      histsVBFPreselDiMuPtL20.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      histsVBFPreselDiMuPtL20.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      histsVBFPreselDiMuPtL20.nJets->Fill(mva.nJets, weight);
      histsVBFPreselDiMuPtL20.ht->Fill(mva.ht, weight);
#ifdef BLIND
      if (!(inBlindWindow && isData))
      {
#endif
      histsVBFPreselDiMuPtL20.mDiMu->Fill(mva.mDiMu, weight);
      histsVBFPreselDiMuPtL20.mDiMuResSigUp->Fill(mDiMuResSigUp, weight);
      histsVBFPreselDiMuPtL20.mDiMuResSigDown->Fill(mDiMuResSigDown, weight);
      histsVBFPreselDiMuPtL20.mDiMuResASigUp->Fill(mDiMuResASigUp, weight);
      histsVBFPreselDiMuPtL20.mDiMuResASigDown->Fill(mDiMuResASigDown, weight);

      if(!vbfPreselection)
      {
        histsVBFPreselDiMuPtL20.BDTHistMuonOnly->Fill(bdtValInc, weight);
        histsVBFPreselDiMuPtL20.likelihoodHistMuonOnly->Fill(likeValInc, weight);
        histsVBFPreselDiMuPtL20.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
      }
      else
      {
        histsVBFPreselDiMuPtL20.BDTHistVBF->Fill(bdtValVBF, weight);
        histsVBFPreselDiMuPtL20.likelihoodHistVBF->Fill(likeValVBF, weight);
        histsVBFPreselDiMuPtL20.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
      }
#ifdef BLIND
      }
#endif
    }

    timeReading += difftime(timeStopReading,timeStartReading);
    timeProcessing += difftime(timeStartFilling,timeStopReading);
    timeFilling += difftime(time(NULL),timeStartFilling);
  }// end event loop
  time_t timeEndEventLoop = time(NULL);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  hists.Write(outFile,"");
  histsBB.Write(outFile,"BB");
  histsBO.Write(outFile,"BO");
  histsBE.Write(outFile,"BE");
  histsOO.Write(outFile,"OO");
  histsOE.Write(outFile,"OE");
  histsEE.Write(outFile,"EE");
  histsNotBB.Write(outFile,"NotBB");
  hists4GeVWindow.Write(outFile,"4GeVWindow");

  histsVBFPresel.Write(outFile,"VBFPresel");
  histsVBFPreselBB.Write(outFile,"VBFPreselBB");
  histsVBFPreselNotBB.Write(outFile,"VBFPreselNotBB");

  histsIncPresel.Write(outFile,"IncPresel");
  histsIncPreselBB.Write(outFile,"IncPreselBB");
  histsIncPreselBO.Write(outFile,"IncPreselBO");
  histsIncPreselBE.Write(outFile,"IncPreselBE");
  histsIncPreselOO.Write(outFile,"IncPreselOO");
  histsIncPreselOE.Write(outFile,"IncPreselOE");
  histsIncPreselEE.Write(outFile,"IncPreselEE");
  histsIncPreselNotBB.Write(outFile,"IncPreselNotBB");

  histsIncBDTSig80.Write(outFile,"IncBDTSig80");
  histsIncBDTSig80BB.Write(outFile,"IncBDTSig80BB");
  histsIncBDTSig80BO.Write(outFile,"IncBDTSig80BO");
  histsIncBDTSig80BE.Write(outFile,"IncBDTSig80BE");
  histsIncBDTSig80OO.Write(outFile,"IncBDTSig80OO");
  histsIncBDTSig80OE.Write(outFile,"IncBDTSig80OE");
  histsIncBDTSig80EE.Write(outFile,"IncBDTSig80EE");
  histsIncBDTSig80NotBB.Write(outFile,"IncBDTSig80NotBB");

  histsVBFBDTSig80.Write(outFile,"VBFBDTSig80");
  histsVBFBDTSig80BB.Write(outFile,"VBFBDTSig80BB");
  histsVBFBDTSig80.Write(outFile,"VBFBDTSig80NotBB");

  histsVBFPreselDiMuPtL20.Write(outFile,"VBFPreselDiMuPtL20");
  histsIncPreselDiMuPtL20.Write(outFile,"IncPreselDiMuPtL20");

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

  mDiMuResSigUp = new TH1F("mDiMuResSigUp","DiMuon Mass Systematic Shift Up: Sigma",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuResSigUp);
  mDiMuResSigDown = new TH1F("mDiMuResSigDown","DiMuon Mass Systematic Shift Down: Sigma",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuResSigDown);
  mDiMuResASigUp = new TH1F("mDiMuResASigUp","DiMuon Mass Systematic Shift Up: ASigma",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuResASigUp);
  mDiMuResASigDown = new TH1F("mDiMuResASigDown","DiMuon Mass Systematic Shift Down: ASigma",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuResASigDown);

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
HistStruct::Write(TFile* outfile, std::string directory)
{
  if(directory == "")
  {
    outfile->cd();
  }
  else
  {
    TDirectory* dir = outfile->mkdir(directory.c_str());
    dir->cd();
  }

  std::vector<TH1F*>::iterator hist;
  std::vector<TH2F*>::iterator hist2D;
  for(hist = histVec.begin();hist != histVec.end(); hist++)
    (*hist)->Write();
  for(hist2D = histVec2D.begin();hist2D != histVec2D.end(); hist2D++)
    (*hist2D)->Write();

  outfile->cd();
}
