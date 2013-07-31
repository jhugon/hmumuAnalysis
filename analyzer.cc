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
#include "annaCalibCode/SmearingTool2011.h"

#include "src/ScaleFactors_2012.h"
#include "src/ScaleFactors_2011.h"

#ifdef MEKD_STANDALONE
#include "MEKD_Wrapper.h"
#endif

#include "TEfficiency.h"
#include <algorithm>    // std::min

#define PUREWEIGHT
//#define ISMEAR 1
#define ISMEAR 2
//#define EFFTEST

#ifdef ROCHESTER
#include "rochester/rochcor2012.h"
#include "rochester/rochcor.h"
#endif
#ifdef ROCHESTER2012ReReco
#include "rochester/rochcor2012jan22.h"
#include "rochester/rochcor.h"
#endif
#ifdef MUSCLEFIT
#include "MuScleFitCorrector.h"
#endif

// for check doublets in Gen level:
#include <set>
typedef std::pair<float,float> pairOfDouble;

using namespace std;
using namespace boost;

// 7 subcategories in the analysis + 3 useful selections 
//const int NsubCat = 17;
//const int NsubCat = 23;
const int NsubCat = 44;
TString subCategory[NsubCat] = {"IncPreselPtG10BB             ",
                                "IncPreselPtG10BO             ",
                                "IncPreselPtG10BE             ",
                                "IncPreselPtG10OO             ",
                                "IncPreselPtG10OE             ",
                                "IncPreselPtG10EE             ",
                                "IncPresel                    ",
                                "IncPreselPtG10               ",
                                "VBFPresel                    ",
                                "VBFBDTCut                    ",
                                "IncPreselBB                  ",
                                "IncPreselBO                  ",
                                "IncPreselBE                  ",
                                "IncPreselOO                  ",
                                "IncPreselOE                  ",
                                "IncPreselEE                  ",
                                "VBFMJJG550                   ", 
                                "IncPreselPtG10BBres1         ",
                                "IncPreselPtG10BBres2         ",
                                "IncPreselPtG10BOres1         ",
                                "IncPreselPtG10BOres2         ",
                                "VBFDeJJG3p5MJJG550pTmissL100 ",
                                "VBFDeJJG3p4MJJG500pTmissL25  ",
                                // baseline++ from i = 23
                                "Jets01PassPtG10    ",
                                "Jets01PassPtG10BB  ",
                                "Jets01PassPtG10BO  ",
                                "Jets01PassPtG10BE  ",
                                "Jets01PassPtG10OO  ",
                                "Jets01PassPtG10OE  ",
                                "Jets01PassPtG10EE  ",
                                "Jets01PassPtG10CC  ", // BE+OO
                                "Jets01PassPtG10FF  ", // OE+EE

                                "Jets01FailPtG10    ",
                                "Jets01FailPtG10BB  ",
                                "Jets01FailPtG10BO  ",
                                "Jets01FailPtG10BE  ",
                                "Jets01FailPtG10OO  ",
                                "Jets01FailPtG10OE  ",
                                "Jets01FailPtG10EE  ",
                                "Jets01FailPtG10CC  ", // BE+OO
                                "Jets01FailPtG10FF  ", // OE+EE

                                "Jet2CutsVBFPass    ",
                                "Jet2CutsGFPass     ",
                                "Jet2CutsFailVBFGF  "};
int Flag_subCat = 1; // 1 make sub category eff and don't split on Nsample <- used
// 0 don't make sub category and split on Nsample <- not used now
const int Nsample = 1; // devide event on subsamples <- not used now 

Double_t MASS_MUON = 0.105658367;    //GeV/c2



void fillMuonHist(TH1F* hist, _MuonInfo& mu1, _MuonInfo& mu2);
void printStationMiss(_MuonInfo& mu1, _MuonInfo& mu2, 
                      _EventInfo& eventInfo, std::string & testString, 
                      unsigned & testCounter);
double GetMassRes(_MuonInfo& mu1, _MuonInfo& mu2);

struct HistStruct
{
  HistStruct();
  ~HistStruct();
  void Write(TFile* outfile, std::string directory);
  void Fill(const MVA& mva, bool blind);

  std::vector<TH1F*> histVec;
  std::vector<TH2F*> histVec2D;
 
  TH1F* mDiMu; // default is PFIso
  TH1F* mDiMu110to160; 
  TH1F* mDiMuTrkRelIso;

  TH1F* RelMassRes;
  TH1F* RelMassResCov;

  TH1F* mDiJet;
  TH1F* ptDiMu;
  TH1F* ptDiJet;
  TH1F* yDiMu;
  TH1F* yDiJet;

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
  TH1F* deltaPhiHJ1;

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

  TH1F* BDTHistVBF;

  TH2F* BDTHistMuonOnlyVMass;
  TH2F* BDTHistVBFVMass;

  TH1F* relIsoMu1;
  TH1F* relIsoMu2;

  TH1F* trkRelIsoMu1;
  TH1F* trkRelIsoMu2;

  TH1F* nJets;
  TH1F* ht;
  TH1F* nJetsInRapidityGap;
  TH1F* htInRapidityGap;

  TH1F* nPU;
  TH1F* nVtx;
  TH1F* met;
  TH1F* ptmiss;
  TH1F* weight;

};

int main(int argc, char *argv[])
{
  float minMmm = 70.0;  // 110;//
  float maxMmm = 400.0; // 160;//
  static const float dataMCMinMass = 110.;
  static const float dataMCMaxMass = 150.;
  static const float MinMassEff = 110.;
  static const float MaxMassEff = 160.;
  //static const float MaxMassEff = 150.;
  float minBlind = 120;
  float maxBlind = 130;

  gErrorIgnoreLevel = kError;
  //gErrorIgnoreLevel = kUnset;
  //gDebug = kLogCrit;
  time_t timeStart = time(NULL);


  set<pairOfDouble> uniqueGeneratedEvents;

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
    ("runPeriod,r",program_options::value<string>(), "Running Periods e.g. 7TeV, 8TeV")
    ("dataMC,d", (std::string("Flag to only keep events in mass range: ").appendAny(dataMCMinMass).append(" to ").appendAny(dataMCMaxMass)).c_str())
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

  if (optionMap.count("dataMC"))
    {
      cout << "Only Running in Data/MC Mass range\n";
      minMmm = dataMCMinMass;
      maxMmm = dataMCMaxMass;
    }
  cout << std::string("Cutting on m_{\\mu\\mu} in [").appendAny(minMmm).append(",").appendAny(maxMmm).append("] \n");

  /////////////////////////////
  //////////// Setup //////////
  /////////////////////////////

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
  signalWords.push_back("Hmm");
  signalWords.push_back("HToMM");
  signalWords.push_back("HToMuMu");
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

  // Run periods
  bool is2011A = false;
  bool is2011B = false;
  bool is2012A = false;
  bool is2012B = false;
  bool is2012C = false;
  bool is2012D = false;
  static const regex re2011A("2011A");
  static const regex re2011B("2011B");
  static const regex re2012A("2012A");
  static const regex re2012B("2012B");
  static const regex re2012C("2012C");
  static const regex re2012D("2012D");
  if(regex_search(filenames[0],re2011A))
    is2011A = true;
  if(regex_search(filenames[0],re2011B))
    is2011B = true;
  if(regex_search(filenames[0],re2012A))
    is2012A = true;
  if(regex_search(filenames[0],re2012B))
    is2012B = true;
  if(regex_search(filenames[0],re2012C))
    is2012C = true;
  if(regex_search(filenames[0],re2012D))
    is2012D = true;
  if(is2011A)
    std::cout << "This is 2011A\n";
  if(is2011B)
    std::cout << "This is 2011B\n";
  if(is2012A)
    std::cout << "This is 2012A\n";
  if(is2012B)
    std::cout << "This is 2012B\n";
  if(is2012C)
    std::cout << "This is 2012C\n";
  if(is2012D)
    std::cout << "This is 2012D\n";

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
      muonIdFuncPtr = &isKinTight_2011_noIso;
    }
  else
    {
      cout << "Using 2012 Tight Muon Selection\n";
      muonIdFuncPtr = &isKinTight_2012_noIso;
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
  float recoCandMassRes, recoCandMassResCov;

  tree->SetBranchAddress("recoCandMass",       &recoCandMass);
  tree->SetBranchAddress("recoCandPt",         &recoCandPt);
  tree->SetBranchAddress("recoCandY",          &recoCandY);
  tree->SetBranchAddress("recoCandPhi",        &recoCandPhi);
  tree->SetBranchAddress("recoCandMassRes",    &recoCandMassRes);
  tree->SetBranchAddress("recoCandMassResCov", &recoCandMassResCov);

  float trueMass=-99999.0;
  if(!isData && tree->GetBranchStatus("trueMass"))
    tree->SetBranchAddress("trueMass", &trueMass);

  /// Higgs Boson 
  _genPartInfo genHpostFSR;
  if(!isData && tree->GetBranchStatus("genHpostFSR"))
    tree->SetBranchAddress("genHpostFSR", &genHpostFSR);

  _TrackInfo reco1GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM1HpostFSR"))
    tree->SetBranchAddress("genM1HpostFSR", &reco1GenPostFSR);

  _TrackInfo reco2GenPostFSR;
  if(!isData && tree->GetBranchStatus("genM2HpostFSR"))
    tree->SetBranchAddress("genM2HpostFSR", &reco2GenPostFSR);

  _PFJetInfo jets;
  tree->SetBranchAddress("pfJets",&jets);

  float puJetFullDisc[10];
  float puJetSimpleDisc[10];
  float puJetCutDisc[10];

  tree->SetBranchAddress("puJetFullDisc",&puJetFullDisc);
  tree->SetBranchAddress("puJetSimpleDisc",&puJetSimpleDisc);
  tree->SetBranchAddress("puJetCutDisc",&puJetCutDisc);

  float puJetFullId[10];
  float puJetSimpleId[10];
  float puJetCutId[10];

  tree->SetBranchAddress("puJetFullId",&puJetFullId);
  tree->SetBranchAddress("puJetSimpleId",&puJetSimpleId);
  tree->SetBranchAddress("puJetCutId",&puJetCutId);

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

  HistStruct histsIncPresel; //
  HistStruct histsIncPreselBB;
  HistStruct histsIncPreselBO;
  HistStruct histsIncPreselBE;
  HistStruct histsIncPreselOO;
  HistStruct histsIncPreselOE;
  HistStruct histsIncPreselEE;
  HistStruct histsIncPreselNotBB;

  HistStruct histsVBFPresel; //
  HistStruct histsVBFPreselBB;
  HistStruct histsVBFPreselNotBB;

  HistStruct histsIncBDTCut;      //     
  HistStruct histsIncBDTCutBB;    //
  HistStruct histsIncBDTCutBO;    //
  HistStruct histsIncBDTCutBE;    //
  HistStruct histsIncBDTCutOO;    //
  HistStruct histsIncBDTCutOE;    //
  HistStruct histsIncBDTCutEE;    //
  HistStruct histsIncBDTCutNotBB; //
                                  
  HistStruct histsVBFBDTCut;      //
  HistStruct histsVBFBDTCutBB;    //   
  HistStruct histsVBFBDTCutNotBB; //

  HistStruct histsIncPreselDiMuPtL20;
  HistStruct histsVBFPreselDiMuPtL20;

  HistStruct histsIncPreselPUJETID;
  HistStruct histsVBFPreselPUJETID;

  HistStruct histsIncPreselPUJETIDForVeto;
  HistStruct histsVBFPreselPUJETIDForVeto;

  HistStruct histsVBFPreselPtMiss50Veto;

  HistStruct histsIncPreselPtG10;
  HistStruct histsIncPreselPtG10BB;
  HistStruct histsIncPreselPtG10BO;
  HistStruct histsIncPreselPtG10BE;
  HistStruct histsIncPreselPtG10OO;
  HistStruct histsIncPreselPtG10OE;
  HistStruct histsIncPreselPtG10EE;
  HistStruct histsIncPreselPtG10NotBB;
  HistStruct histsIncPreselPtG10BBres1;
  HistStruct histsIncPreselPtG10BBres2;
  HistStruct histsIncPreselPtG10BOres1;
  HistStruct histsIncPreselPtG10BOres2;
  HistStruct histsIncPreselPtG10BBres1Cov;
  HistStruct histsIncPreselPtG10BBres2Cov;

  HistStruct histsVBFMJJG550;
  // couple of VBF CiC optimization points
  HistStruct histsVBFDeJJG3p5MJJG550pTmissL100;
  HistStruct histsVBFDeJJG3p4MJJG500pTmissL25;

  // baseline++
  HistStruct histsJets01PassPtG10;
  HistStruct histsJets01PassPtG10BB;
  HistStruct histsJets01PassPtG10BO;
  HistStruct histsJets01PassPtG10BE;
  HistStruct histsJets01PassPtG10OO;
  HistStruct histsJets01PassPtG10OE;
  HistStruct histsJets01PassPtG10EE;
  HistStruct histsJets01PassPtG10CC; // BE+OO
  HistStruct histsJets01PassPtG10FF; // OE+EE

  HistStruct histsJets01FailPtG10;
  HistStruct histsJets01FailPtG10BB;
  HistStruct histsJets01FailPtG10BO;
  HistStruct histsJets01FailPtG10BE;
  HistStruct histsJets01FailPtG10OO;
  HistStruct histsJets01FailPtG10OE;
  HistStruct histsJets01FailPtG10EE;
  HistStruct histsJets01FailPtG10CC; // BE+OO
  HistStruct histsJets01FailPtG10FF; // OE+EE

  HistStruct histsJet2CutsVBFPass;
  HistStruct histsJet2CutsGFPass;
  HistStruct histsJet2CutsFailVBFGF;

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
  reweight::LumiReWeighting lumiWeights("pileupDists/PileUpHistMC2012Summer50nsPoissonOOTPU.root","pileupDists/PileUpHist2012ABCD.root","pileup","pileup");
  if (runPeriod == "7TeV")
    {
      cout << "Using 2011AB PU reweighting\n";
      lumiWeights = reweight::LumiReWeighting("pileupDists/PileUpHistMCFall11.root","pileupDists/PileUpHist2011AB.root","pileup","pileup");
    }
  else
    {
      cout << "Using 2012ABCD PU reweighting\n";
    }
#endif

  /////////////////////////////
  // Muon Momentum Corrections

#ifdef MUSCLEFIT
  MuScleFitCorrector* mfCorr;
  TString mfInFile;
#ifdef MUSCLEFIT2012ReReco
  if (isData)
  {
    if(runPeriod == "8TeV")
    {
      if(is2012D)
        mfInFile = "musclefitReReco/MuScleFit_2012D_DATA_ReReco_53X.txt";
      else
        mfInFile = "musclefitReReco/MuScleFit_2012ABC_DATA_ReReco_53X.txt";
    }
    else
      mfInFile = "musclefitReReco/MuScleFit_2011_DATA_44X.txt";
  }
  else
  {
    if(runPeriod == "8TeV")
      mfInFile = "musclefitReReco/MuScleFit_2012_MC_53X_smearReReco.txt";
    else
      mfInFile = "musclefitReReco/MuScleFit_2011_MC_44X.txt";
  }
#else
  if (isData)
    {
      if(runPeriod == "8TeV")
        mfInFile = "musclefit/MuScleFit_2012_DATA_53X.txt";
      else
        mfInFile = "musclefit/MuScleFit_2011_DATA_44X.txt";
    }
  else
  {
    if(runPeriod == "8TeV")
      mfInFile = "musclefit/MuScleFit_2012_MC_52X.txt";
    else
      mfInFile = "musclefit/MuScleFit_2011_MC_44X.txt";
  }
#endif
  cout << "Using musclefit file: "<<mfInFile<<endl;
  mfCorr = new MuScleFitCorrector(mfInFile);
#endif
#ifdef ROCHESTER
  rochcor2012* rCorr12 = new rochcor2012();
  rochcor* rCorr11 = new rochcor();
  int rochesterRun=0;
  if(is2011B)
    rochesterRun=1;
#endif
#ifdef ROCHESTER2012ReReco
  rochcor2012jan22* rCorr12jan22 = new rochcor2012jan22();
#endif

  /////////////////////////
  // Smearing
  SmearingTool *smearPT = new SmearingTool();
  SmearingTool2011 *smearPT2011 = new SmearingTool2011();

  const double SQRT2 = sqrt(2);

  TRandom3 randomForSF(123412845);

  //////////////////////////
  // Defining the outTree
  outFile->cd();
  TTree* _outTree = new TTree("outtree", "myTree");
  cout << "outtree directory name: " << _outTree->GetDirectory()->GetName()<<endl;

  int _eventInfo_run;
  int _eventInfo_lumi;
  int _eventInfo_event;

  float _nPU;
  float _puWeight;
  float _puWeightMuonEffUp;
  float _puWeightMuonEffDown;
  int _eventType;

  float _dimuonMass; 
  float _dimuonMassResSFUp; 
  float _dimuonMassResSFDown; 
  float _dimuonPt; 
  float _dimuonY;    
  float _cosThetaStar;

  int   _muonLead_charge;      
  float _muonLead_pt;          
  float _muonLead_ptErr;       
  float _muonLead_eta;         
  float _muonLead_phi;         
  int   _muonLead_passTrkRelIso;
  int   _muonLead_passPFRelIso;
  int   _muonLead_isHltMatched;

  int   _muonSub_charge;      
  float _muonSub_pt;          
  float _muonSub_ptErr;       
  float _muonSub_eta;         
  float _muonSub_phi;         
  int   _muonSub_passTrkRelIso;
  int   _muonSub_passPFRelIso;
  int   _muonSub_isHltMatched;
                                                                    
  float _nJets;
  float _ptMiss;
  float _deltaEtaJets;
  float _dijetMass;
  float _dijetPt;
  float _dijetY;

  float _jetLead_pt;          
  float _jetLead_eta;         
  float _jetLead_PUIDDisc;         
  int   _jetLead_PUIDFlag;         
  float _jetLead_FullPUIDDisc;         
  int   _jetLead_FullPUIDFlag;         
  float _jetLead_CutPUIDDisc;         
  int   _jetLead_flav;          

  float _jetSub_pt;          
  float _jetSub_eta;         
  float _jetSub_PUIDDisc;         
  int   _jetSub_PUIDFlag;         
  float _jetSub_FullPUIDDisc;         
  int   _jetSub_FullPUIDFlag;         
  float _jetSub_CutPUIDDisc;         
  int   _jetSub_flav;          

  float _kd;
  float _sigME;
  float _bakME;
  float _kdPdf;
  float _sigMEPdf;
  float _bakMEPdf;
  float _bdtNonVBF;
  float _bdtVBF;
     
  // eventInfo;
  _outTree->Branch("eventInfo_run",    &_eventInfo_run,   "eventInfo_run/I");
  _outTree->Branch("eventInfo_lumi",   &_eventInfo_lumi,  "eventInfo_lumi/I");
  _outTree->Branch("eventInfo_event",  &_eventInfo_event, "eventInfo_event/I");

  _outTree->Branch("nPU",                 &_nPU,                 "nPU/F");
  _outTree->Branch("puWeight",            &_puWeight,            "puWeight/F");
  _outTree->Branch("puWeightMuonEffUp",   &_puWeightMuonEffUp,   "puWeightMuonEffUp/F");
  _outTree->Branch("puWeightMuonEffDown", &_puWeightMuonEffDown, "puWeightMuonEffDown/F");

  _outTree->Branch("eventType",    &_eventType,   "eventType/I");

  _outTree->Branch("dimuonMass",            &_dimuonMass,          "dimuonMass/F");
  _outTree->Branch("dimuonMassResSFUp",     &_dimuonMassResSFUp,   "dimuonMassResSFUp/F");
  _outTree->Branch("dimuonMassResSFDown",   &_dimuonMassResSFDown, "dimuonMassResSFDown/F");
  _outTree->Branch("dimuonPt"  ,   &_dimuonPt  ,   "dimuonPt/F");
  _outTree->Branch("dimuonY",      &_dimuonY,      "dimuonY/F");
  _outTree->Branch("cosThetaStar", &_cosThetaStar, "cosThetaStar/F");

  _outTree->Branch("muonLead_charge",       &_muonLead_charge,         "muonLead_charge/I");
  _outTree->Branch("muonLead_pt",           &_muonLead_pt,             "muonLead_pt/F");          
  _outTree->Branch("muonLead_ptErr",        &_muonLead_ptErr,          "muonLead_ptErr/F");       
  _outTree->Branch("muonLead_eta",          &_muonLead_eta,            "muonLead_eta/F");         
  _outTree->Branch("muonLead_phi",          &_muonLead_phi,            "muonLead_phi/F");         
  _outTree->Branch("muonLead_passTrkRelIso",&_muonLead_passTrkRelIso,  "muonLead_passTrkRelIso/I"); 
  _outTree->Branch("muonLead_passPFRelIso", &_muonLead_passPFRelIso,   "muonLead_passPFRelIso/I");
  _outTree->Branch("muonLead_isHltMatched", &_muonLead_isHltMatched,   "muonLead_isHltMatched/I");

  _outTree->Branch("muonSub_charge",       &_muonSub_charge,         "muonSub_charge/I");
  _outTree->Branch("muonSub_pt",           &_muonSub_pt,             "muonSub_pt/F");          
  _outTree->Branch("muonSub_ptErr",        &_muonSub_ptErr,          "muonSub_ptErr/F");       
  _outTree->Branch("muonSub_eta",          &_muonSub_eta,            "muonSub_eta/F");         
  _outTree->Branch("muonSub_phi",          &_muonSub_phi,            "muonSub_phi/F");         
  _outTree->Branch("muonSub_passTrkRelIso",&_muonSub_passTrkRelIso,  "muonSub_passTrkRelIso/I"); 
  _outTree->Branch("muonSub_passPFRelIso", &_muonSub_passPFRelIso,   "muonSub_passPFRelIso/I");
  _outTree->Branch("muonSub_isHltMatched", &_muonSub_isHltMatched,   "muonSub_isHltMatched/I");

  _outTree->Branch("nJets",        &_nJets,        "nJets/F");
  _outTree->Branch("ptMiss",       &_ptMiss,       "ptMiss/F");
  _outTree->Branch("deltaEtaJets", &_deltaEtaJets, "deltaEtaJets/F");
  _outTree->Branch("dijetMass",       &_dijetMass,       "dijetMass/F");
  _outTree->Branch("dijetPt",       &_dijetPt,       "dijetPt/F");
  _outTree->Branch("dijetY",       &_dijetY,       "dijetY/F");

  _outTree->Branch("jetLead_pt",           &_jetLead_pt,             "jetLead_pt/F");          
  _outTree->Branch("jetLead_eta",          &_jetLead_eta,            "jetLead_eta/F");         
  _outTree->Branch("jetLead_PUIDDisc",          &_jetLead_PUIDDisc,            "jetLead_PUIDDisc/F");         
  _outTree->Branch("jetLead_PUIDFlag",          &_jetLead_PUIDFlag,            "jetLead_PUIDFlag/I");         
  _outTree->Branch("jetLead_FullPUIDDisc",          &_jetLead_FullPUIDDisc,            "jetLead_FullPUIDDisc/F");         
  _outTree->Branch("jetLead_FullPUIDFlag",          &_jetLead_FullPUIDFlag,            "jetLead_FullPUIDFlag/I");         
  _outTree->Branch("jetLead_CutPUIDDisc",          &_jetLead_CutPUIDDisc,            "jetLead_CutPUIDDisc/F");         
  _outTree->Branch("jetLead_flav",           &_jetLead_flav,             "jetLead_flav/I");          

  _outTree->Branch("jetSub_pt",           &_jetSub_pt,             "jetSub_pt/F");          
  _outTree->Branch("jetSub_eta",          &_jetSub_eta,            "jetSub_eta/F");         
  _outTree->Branch("jetSub_PUIDDisc",          &_jetSub_PUIDDisc,            "jetSub_PUIDDisc/F");         
  _outTree->Branch("jetSub_PUIDFlag",          &_jetSub_PUIDFlag,            "jetSub_PUIDFlag/I");         
  _outTree->Branch("jetSub_FullPUIDDisc",          &_jetSub_FullPUIDDisc,            "jetSub_FullPUIDDisc/F");         
  _outTree->Branch("jetSub_FullPUIDFlag",          &_jetSub_FullPUIDFlag,            "jetSub_FullPUIDFlag/I");         
  _outTree->Branch("jetSub_CutPUIDDisc",          &_jetSub_CutPUIDDisc,            "jetSub_CutPUIDDisc/F");         
  _outTree->Branch("jetSub_flav",           &_jetSub_flav,             "jetSub_flav/I");          

  _outTree->Branch("kd",           &_kd,           "kd/F");
  _outTree->Branch("sigME",        &_sigME,        "sigME/F");
  _outTree->Branch("bakME",        &_bakME,        "bakME/F");
  _outTree->Branch("kdPdf",           &_kdPdf,           "kdPdf/F");
  _outTree->Branch("sigMEPdf",        &_sigMEPdf,        "sigMEPdf/F");
  _outTree->Branch("bakMEPdf",        &_bakMEPdf,        "bakMEPdf/F");
  _outTree->Branch("bdtVBF",        &_bdtVBF,        "bdtVBF/F");
  _outTree->Branch("bdtNonVBF",        &_bdtNonVBF,        "bdtNonVBF/F");

  float _dimuonMass_noMuscle;
  float _dimuonPt_noMuscle;
  float _muonLead_pt_noMuscle;
  float _muonSub_pt_noMuscle;
  int   _muonLead_passTrkRelIso_noMuscle;
  int   _muonLead_passPFRelIso_noMuscle;
  int   _muonSub_passTrkRelIso_noMuscle;
  int   _muonSub_passPFRelIso_noMuscle;

  _outTree->Branch("dimuonMass_noMuscle",   &_dimuonMass_noMuscle,   "dimuonMass_noMuscle/F");
  _outTree->Branch("dimuonPt_noMuscle"  ,   &_dimuonPt_noMuscle  ,   "dimuonPt_noMuscle/F");
  _outTree->Branch("muonLead_pt_noMuscle",           &_muonLead_pt_noMuscle,             "muonLead_pt_noMuscle/F");          
  _outTree->Branch("muonSub_pt_noMuscle",           &_muonSub_pt_noMuscle,             "muonSub_pt_noMuscle/F");          
  _outTree->Branch("muonLead_passTrkRelIso_noMuscle",&_muonLead_passTrkRelIso_noMuscle,  "muonLead_passTrkRelIso_noMuscle/I"); 
  _outTree->Branch("muonLead_passPFRelIso_noMuscle", &_muonLead_passPFRelIso_noMuscle,   "muonLead_passPFRelIso_noMuscle/I");
  _outTree->Branch("muonSub_passTrkRelIso_noMuscle",&_muonSub_passTrkRelIso_noMuscle,  "muonSub_passTrkRelIso_noMuscle/I"); 
  _outTree->Branch("muonSub_passPFRelIso_noMuscle", &_muonSub_passPFRelIso_noMuscle,   "muonSub_passPFRelIso_noMuscle/I");

  float _nJets_JESUp;
  float _ptMiss_JESUp;
  float _deltaEtaJets_JESUp;
  float _dijetMass_JESUp;

  float _jetLead_pt_JESUp;          
  float _jetLead_eta_JESUp;         
  float _jetSub_pt_JESUp;          
  float _jetSub_eta_JESUp;         

  _outTree->Branch("nJets_JESUp",        &_nJets_JESUp,        "nJets_JESUp/F");
  _outTree->Branch("ptMiss_JESUp",       &_ptMiss_JESUp,       "ptMiss_JESUp/F");
  _outTree->Branch("deltaEtaJets_JESUp", &_deltaEtaJets_JESUp, "deltaEtaJets_JESUp/F");
  _outTree->Branch("dijetMass_JESUp",       &_dijetMass_JESUp,       "dijetMass_JESUp/F");
  _outTree->Branch("jetLead_pt_JESUp",           &_jetLead_pt_JESUp,             "jetLead_pt_JESUp/F");          
  _outTree->Branch("jetLead_eta_JESUp",          &_jetLead_eta_JESUp,            "jetLead_eta_JESUp/F");         
  _outTree->Branch("jetSub_pt_JESUp",           &_jetSub_pt_JESUp,             "jetSub_pt_JESUp/F");          
  _outTree->Branch("jetSub_eta_JESUp",          &_jetSub_eta_JESUp,            "jetSub_eta_JESUp/F");         

  float _nJets_JESDown;
  float _ptMiss_JESDown;
  float _deltaEtaJets_JESDown;
  float _dijetMass_JESDown;

  float _jetLead_pt_JESDown;          
  float _jetLead_eta_JESDown;         
  float _jetSub_pt_JESDown;          
  float _jetSub_eta_JESDown;         

  _outTree->Branch("nJets_JESDown",        &_nJets_JESDown,        "nJets_JESDown/F");
  _outTree->Branch("ptMiss_JESDown",       &_ptMiss_JESDown,       "ptMiss_JESDown/F");
  _outTree->Branch("deltaEtaJets_JESDown", &_deltaEtaJets_JESDown, "deltaEtaJets_JESDown/F");
  _outTree->Branch("dijetMass_JESDown",       &_dijetMass_JESDown,       "dijetMass_JESDown/F");
  _outTree->Branch("jetLead_pt_JESDown",           &_jetLead_pt_JESDown,             "jetLead_pt_JESDown/F");          
  _outTree->Branch("jetLead_eta_JESDown",          &_jetLead_eta_JESDown,            "jetLead_eta_JESDown/F");         
  _outTree->Branch("jetSub_pt_JESDown",           &_jetSub_pt_JESDown,             "jetSub_pt_JESDown/F");          
  _outTree->Branch("jetSub_eta_JESDown",          &_jetSub_eta_JESDown,            "jetSub_eta_JESDown/F");         

  float _nJets_JERUp;
  float _ptMiss_JERUp;
  float _deltaEtaJets_JERUp;
  float _dijetMass_JERUp;

  float _jetLead_pt_JERUp;          
  float _jetLead_eta_JERUp;         
  float _jetSub_pt_JERUp;          
  float _jetSub_eta_JERUp;         

  _outTree->Branch("nJets_JERUp",        &_nJets_JERUp,        "nJets_JERUp/F");
  _outTree->Branch("ptMiss_JERUp",       &_ptMiss_JERUp,       "ptMiss_JERUp/F");
  _outTree->Branch("deltaEtaJets_JERUp", &_deltaEtaJets_JERUp, "deltaEtaJets_JERUp/F");
  _outTree->Branch("dijetMass_JERUp",       &_dijetMass_JERUp,       "dijetMass_JERUp/F");
  _outTree->Branch("jetLead_pt_JERUp",           &_jetLead_pt_JERUp,             "jetLead_pt_JERUp/F");          
  _outTree->Branch("jetLead_eta_JERUp",          &_jetLead_eta_JERUp,            "jetLead_eta_JERUp/F");         
  _outTree->Branch("jetSub_pt_JERUp",           &_jetSub_pt_JERUp,             "jetSub_pt_JERUp/F");          
  _outTree->Branch("jetSub_eta_JERUp",          &_jetSub_eta_JERUp,            "jetSub_eta_JERUp/F");         

  float _nJets_JERDown;
  float _ptMiss_JERDown;
  float _deltaEtaJets_JERDown;
  float _dijetMass_JERDown;

  float _jetLead_pt_JERDown;          
  float _jetLead_eta_JERDown;         
  float _jetSub_pt_JERDown;          
  float _jetSub_eta_JERDown;         

  _outTree->Branch("nJets_JERDown",        &_nJets_JERDown,        "nJets_JERDown/F");
  _outTree->Branch("ptMiss_JERDown",       &_ptMiss_JERDown,       "ptMiss_JERDown/F");
  _outTree->Branch("deltaEtaJets_JERDown", &_deltaEtaJets_JERDown, "deltaEtaJets_JERDown/F");
  _outTree->Branch("dijetMass_JERDown",       &_dijetMass_JERDown,       "dijetMass_JERDown/F");
  _outTree->Branch("jetLead_pt_JERDown",           &_jetLead_pt_JERDown,             "jetLead_pt_JERDown/F");          
  _outTree->Branch("jetLead_eta_JERDown",          &_jetLead_eta_JERDown,            "jetLead_eta_JERDown/F");         
  _outTree->Branch("jetSub_pt_JERDown",           &_jetSub_pt_JERDown,             "jetSub_pt_JERDown/F");          
  _outTree->Branch("jetSub_eta_JERDown",          &_jetSub_eta_JERDown,            "jetSub_eta_JERDown/F");         

  float _puidUncWeight;         
  _outTree->Branch("puidUncWeight",     &_puidUncWeight,    "puidUncWeight/F");

  // Defining the outTree
  TTree* treeJets01PassPtG10   = new TTree("outtreeJets01PassPtG10","");
  TTree* treeJets01PassPtG10BB = new TTree("outtreeJets01PassPtG10BB","");
  TTree* treeJets01PassPtG10BO = new TTree("outtreeJets01PassPtG10BO","");
  TTree* treeJets01PassPtG10BE = new TTree("outtreeJets01PassPtG10BE","");
  TTree* treeJets01PassPtG10OO = new TTree("outtreeJets01PassPtG10OO","");
  TTree* treeJets01PassPtG10OE = new TTree("outtreeJets01PassPtG10OE","");
  TTree* treeJets01PassPtG10EE = new TTree("outtreeJets01PassPtG10EE","");
                                            
  TTree* treeJets01FailPtG10   = new TTree("outtreeJets01FailPtG10","");
  TTree* treeJets01FailPtG10BB = new TTree("outtreeJets01FailPtG10BB","");
  TTree* treeJets01FailPtG10BO = new TTree("outtreeJets01FailPtG10BO","");
  TTree* treeJets01FailPtG10BE = new TTree("outtreeJets01FailPtG10BE","");
  TTree* treeJets01FailPtG10OO = new TTree("outtreeJets01FailPtG10OO","");
  TTree* treeJets01FailPtG10OE = new TTree("outtreeJets01FailPtG10OE","");
  TTree* treeJets01FailPtG10EE = new TTree("outtreeJets01FailPtG10EE","");
                                            
  TTree* treeJet2CutsVBFPass   = new TTree("outtreeJet2CutsVBFPass","");
  TTree* treeJet2CutsGFPass    = new TTree("outtreeJet2CutsGFPass","");
  TTree* treeJet2CutsFailVBFGF = new TTree("outtreeJet2CutsFailVBFGF","");


  treeJets01PassPtG10   -> SetDirectory(0);
  treeJets01PassPtG10BB -> SetDirectory(0);
  treeJets01PassPtG10BO -> SetDirectory(0);
  treeJets01PassPtG10BE -> SetDirectory(0);
  treeJets01PassPtG10OO -> SetDirectory(0);
  treeJets01PassPtG10OE -> SetDirectory(0);
  treeJets01PassPtG10EE -> SetDirectory(0);
                      
  treeJets01FailPtG10   -> SetDirectory(0);
  treeJets01FailPtG10BB -> SetDirectory(0);
  treeJets01FailPtG10BO -> SetDirectory(0);
  treeJets01FailPtG10BE -> SetDirectory(0);
  treeJets01FailPtG10OO -> SetDirectory(0);
  treeJets01FailPtG10OE -> SetDirectory(0);
  treeJets01FailPtG10EE -> SetDirectory(0);

  treeJet2CutsVBFPass   -> SetDirectory(0);
  treeJet2CutsGFPass    -> SetDirectory(0);
  treeJet2CutsFailVBFGF -> SetDirectory(0);

  // I am so sad I could not find a better way to do that
  float _dimuonMassJets01PassPtG10  ;  float _dimuonMassJets01PassPtG10ResSFUp  ;  float _dimuonMassJets01PassPtG10ResSFDown  ;
  float _dimuonMassJets01PassPtG10BB;  float _dimuonMassJets01PassPtG10BBResSFUp;  float _dimuonMassJets01PassPtG10BBResSFDown;
  float _dimuonMassJets01PassPtG10BO;  float _dimuonMassJets01PassPtG10BOResSFUp;  float _dimuonMassJets01PassPtG10BOResSFDown;
  float _dimuonMassJets01PassPtG10BE;  float _dimuonMassJets01PassPtG10BEResSFUp;  float _dimuonMassJets01PassPtG10BEResSFDown;
  float _dimuonMassJets01PassPtG10OO;  float _dimuonMassJets01PassPtG10OOResSFUp;  float _dimuonMassJets01PassPtG10OOResSFDown;
  float _dimuonMassJets01PassPtG10OE;  float _dimuonMassJets01PassPtG10OEResSFUp;  float _dimuonMassJets01PassPtG10OEResSFDown;
  float _dimuonMassJets01PassPtG10EE;  float _dimuonMassJets01PassPtG10EEResSFUp;  float _dimuonMassJets01PassPtG10EEResSFDown;

  float _dimuonMassJets01FailPtG10  ;  float _dimuonMassJets01FailPtG10ResSFUp  ;  float _dimuonMassJets01FailPtG10ResSFDown  ;
  float _dimuonMassJets01FailPtG10BB;  float _dimuonMassJets01FailPtG10BBResSFUp;  float _dimuonMassJets01FailPtG10BBResSFDown;
  float _dimuonMassJets01FailPtG10BO;  float _dimuonMassJets01FailPtG10BOResSFUp;  float _dimuonMassJets01FailPtG10BOResSFDown;
  float _dimuonMassJets01FailPtG10BE;  float _dimuonMassJets01FailPtG10BEResSFUp;  float _dimuonMassJets01FailPtG10BEResSFDown;
  float _dimuonMassJets01FailPtG10OO;  float _dimuonMassJets01FailPtG10OOResSFUp;  float _dimuonMassJets01FailPtG10OOResSFDown;
  float _dimuonMassJets01FailPtG10OE;  float _dimuonMassJets01FailPtG10OEResSFUp;  float _dimuonMassJets01FailPtG10OEResSFDown;
  float _dimuonMassJets01FailPtG10EE;  float _dimuonMassJets01FailPtG10EEResSFUp;  float _dimuonMassJets01FailPtG10EEResSFDown;

  float _dimuonMassJet2CutsVBFPass  ;  float _dimuonMassJet2CutsVBFPassResSFUp  ;  float _dimuonMassJet2CutsVBFPassResSFDown  ;
  float _dimuonMassJet2CutsGFPass   ;  float _dimuonMassJet2CutsGFPassResSFUp   ;  float _dimuonMassJet2CutsGFPassResSFDown   ;
  float _dimuonMassJet2CutsFailVBFGF;  float _dimuonMassJet2CutsFailVBFGFResSFUp;  float _dimuonMassJet2CutsFailVBFGFResSFDown;


  treeJets01PassPtG10   -> Branch("dimuonMass",   &_dimuonMassJets01PassPtG10  ,   "dimuonMass/F");
  treeJets01PassPtG10BB -> Branch("dimuonMass",   &_dimuonMassJets01PassPtG10BB,   "dimuonMass/F");
  treeJets01PassPtG10BO -> Branch("dimuonMass",   &_dimuonMassJets01PassPtG10BO,   "dimuonMass/F");
  treeJets01PassPtG10BE -> Branch("dimuonMass",   &_dimuonMassJets01PassPtG10BE,   "dimuonMass/F");
  treeJets01PassPtG10OO -> Branch("dimuonMass",   &_dimuonMassJets01PassPtG10OO,   "dimuonMass/F");
  treeJets01PassPtG10OE -> Branch("dimuonMass",   &_dimuonMassJets01PassPtG10OE,   "dimuonMass/F");
  treeJets01PassPtG10EE -> Branch("dimuonMass",   &_dimuonMassJets01PassPtG10EE,   "dimuonMass/F");
                                                                                     
  treeJets01FailPtG10   -> Branch("dimuonMass",   &_dimuonMassJets01FailPtG10  ,   "dimuonMass/F");
  treeJets01FailPtG10BB -> Branch("dimuonMass",   &_dimuonMassJets01FailPtG10BB,   "dimuonMass/F");
  treeJets01FailPtG10BO -> Branch("dimuonMass",   &_dimuonMassJets01FailPtG10BO,   "dimuonMass/F");
  treeJets01FailPtG10BE -> Branch("dimuonMass",   &_dimuonMassJets01FailPtG10BE,   "dimuonMass/F");
  treeJets01FailPtG10OO -> Branch("dimuonMass",   &_dimuonMassJets01FailPtG10OO,   "dimuonMass/F");
  treeJets01FailPtG10OE -> Branch("dimuonMass",   &_dimuonMassJets01FailPtG10OE,   "dimuonMass/F");
  treeJets01FailPtG10EE -> Branch("dimuonMass",   &_dimuonMassJets01FailPtG10EE,   "dimuonMass/F");
                                                                                     
  treeJet2CutsVBFPass   -> Branch("dimuonMass",   &_dimuonMassJet2CutsVBFPass  ,   "dimuonMass/F");
  treeJet2CutsGFPass    -> Branch("dimuonMass",   &_dimuonMassJet2CutsGFPass   ,   "dimuonMass/F");
  treeJet2CutsFailVBFGF -> Branch("dimuonMass",   &_dimuonMassJet2CutsFailVBFGF,   "dimuonMass/F");


  // Muon Momentum Resolution SF Up
  treeJets01PassPtG10   -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01PassPtG10ResSFUp  ,   "dimuonMassResSFUp/F");
  treeJets01PassPtG10BB -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01PassPtG10BBResSFUp,   "dimuonMassResSFUp/F");
  treeJets01PassPtG10BO -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01PassPtG10BOResSFUp,   "dimuonMassResSFUp/F");
  treeJets01PassPtG10BE -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01PassPtG10BEResSFUp,   "dimuonMassResSFUp/F");
  treeJets01PassPtG10OO -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01PassPtG10OOResSFUp,   "dimuonMassResSFUp/F");
  treeJets01PassPtG10OE -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01PassPtG10OEResSFUp,   "dimuonMassResSFUp/F");
  treeJets01PassPtG10EE -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01PassPtG10EEResSFUp,   "dimuonMassResSFUp/F");

  treeJets01FailPtG10   -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01FailPtG10ResSFUp  ,   "dimuonMassResSFUp/F");
  treeJets01FailPtG10BB -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01FailPtG10BBResSFUp,   "dimuonMassResSFUp/F");
  treeJets01FailPtG10BO -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01FailPtG10BOResSFUp,   "dimuonMassResSFUp/F");
  treeJets01FailPtG10BE -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01FailPtG10BEResSFUp,   "dimuonMassResSFUp/F");
  treeJets01FailPtG10OO -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01FailPtG10OOResSFUp,   "dimuonMassResSFUp/F");
  treeJets01FailPtG10OE -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01FailPtG10OEResSFUp,   "dimuonMassResSFUp/F");
  treeJets01FailPtG10EE -> Branch("dimuonMassResSFUp",   &_dimuonMassJets01FailPtG10EEResSFUp,   "dimuonMassResSFUp/F");

  treeJet2CutsVBFPass   -> Branch("dimuonMassResSFUp",   &_dimuonMassJet2CutsVBFPassResSFUp  ,   "dimuonMassResSFUp/F");
  treeJet2CutsGFPass    -> Branch("dimuonMassResSFUp",   &_dimuonMassJet2CutsGFPassResSFUp   ,   "dimuonMassResSFUp/F");
  treeJet2CutsFailVBFGF -> Branch("dimuonMassResSFUp",   &_dimuonMassJet2CutsFailVBFGFResSFUp,   "dimuonMassResSFUp/F");


  // Muon Momentum Resolution SF Down
  treeJets01PassPtG10   -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01PassPtG10ResSFDown  ,   "dimuonMassResSFDown/F");
  treeJets01PassPtG10BB -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01PassPtG10BBResSFDown,   "dimuonMassResSFDown/F");
  treeJets01PassPtG10BO -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01PassPtG10BOResSFDown,   "dimuonMassResSFDown/F");
  treeJets01PassPtG10BE -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01PassPtG10BEResSFDown,   "dimuonMassResSFDown/F");
  treeJets01PassPtG10OO -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01PassPtG10OOResSFDown,   "dimuonMassResSFDown/F");
  treeJets01PassPtG10OE -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01PassPtG10OEResSFDown,   "dimuonMassResSFDown/F");
  treeJets01PassPtG10EE -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01PassPtG10EEResSFDown,   "dimuonMassResSFDown/F");

  treeJets01FailPtG10   -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01FailPtG10ResSFDown  ,   "dimuonMassResSFDown/F");
  treeJets01FailPtG10BB -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01FailPtG10BBResSFDown,   "dimuonMassResSFDown/F");
  treeJets01FailPtG10BO -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01FailPtG10BOResSFDown,   "dimuonMassResSFDown/F");
  treeJets01FailPtG10BE -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01FailPtG10BEResSFDown,   "dimuonMassResSFDown/F");
  treeJets01FailPtG10OO -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01FailPtG10OOResSFDown,   "dimuonMassResSFDown/F");
  treeJets01FailPtG10OE -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01FailPtG10OEResSFDown,   "dimuonMassResSFDown/F");
  treeJets01FailPtG10EE -> Branch("dimuonMassResSFDown",   &_dimuonMassJets01FailPtG10EEResSFDown,   "dimuonMassResSFDown/F");

  treeJet2CutsVBFPass   -> Branch("dimuonMassResSFDown",   &_dimuonMassJet2CutsVBFPassResSFDown  ,   "dimuonMassResSFDown/F");
  treeJet2CutsGFPass    -> Branch("dimuonMassResSFDown",   &_dimuonMassJet2CutsGFPassResSFDown   ,   "dimuonMassResSFDown/F");
  treeJet2CutsFailVBFGF -> Branch("dimuonMassResSFDown",   &_dimuonMassJet2CutsFailVBFGFResSFDown,   "dimuonMassResSFDown/F");

  //////////////////////////
  // Creating the MEKD

#ifdef MEKD_STANDALONE
  double nTeV = 8.;
  if (runPeriod == "7TeV")
    nTeV = 7.;
  MEKD_Wrapper mekd(nTeV,false);
  MEKD_Wrapper mekdPdf(nTeV,true);
#endif

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  int Nbin;
  if(Flag_subCat == 1) Nbin = NsubCat;
  else Nbin = Nsample;
  //////////////////////////////////////////////////////////////////////
  int counterGenBoson  [Nbin]; // mass cut only
  int counterMinimCuts [Nbin]; // mass, charge, vertex, cosmic
  int counterAccCuts   [Nbin]; // pt, eta
  int counterIDCuts    [Nbin]; // tight muon id
  int counterIsoCuts   [Nbin]; // trk Iso
  int counterTrigCuts  [Nbin]; // iso mu 24_eta2p1
  int counterJetSel    [Nbin]; // jet selection
  int counterBDTvbfCut [Nbin]; // VBF BDT cut
  int counterNonJetSel [Nbin]; // jet selection
  int counterDiPt10GeV [Nbin]; // jet selection

  float EffMinimCuts [Nbin]; // mass, charge, vertex, cosmic
  float EffAccCuts   [Nbin]; // pt, eta
  float EffIDCuts    [Nbin]; // tight muon id
  float EffIsoCuts   [Nbin]; // trk Iso
  float EffTrigCuts  [Nbin]; // iso mu 24_eta2p1
  float EffJetSel    [Nbin]; // jet selection
  float EffBDTvbfCut [Nbin]; // VBF BDT cut
  float EffNonJetSel [Nbin]; // jet selection
  float EffDiPt10GeV [Nbin]; // pt(mumu) > = 10 GeV/c 

  float dEffMinimCuts [Nbin]; // mass, charge, vertex, cosmic
  float dEffAccCuts   [Nbin]; // pt, eta
  float dEffIDCuts    [Nbin]; // tight muon id
  float dEffIsoCuts   [Nbin]; // trk Iso
  float dEffTrigCuts  [Nbin]; // iso mu 24_eta2p1
  float dEffJetSel    [Nbin]; // jet selection
  float dEffBDTvbfCut [Nbin]; // VBF BDT cut
  float dEffNonJetSel [Nbin]; // jet selection
  float dEffDiPt10GeV [Nbin]; // pt(mumu) > = 10 GeV/c 

  //reset them
  for(unsigned iS=0; iS < Nbin;iS++){
    counterGenBoson  [iS] = 0;
    counterMinimCuts [iS] = 0;
    counterAccCuts   [iS] = 0;
    counterIDCuts    [iS] = 0;
    counterIsoCuts   [iS] = 0; // trk Iso
    counterTrigCuts  [iS] = 0; // iso mu 24_eta2p1
    counterJetSel    [iS] = 0; // jet selection
    counterBDTvbfCut [iS] = 0; // VBF BDT cut
    counterNonJetSel [iS] = 0; // jet selection
    counterDiPt10GeV [iS] = 0; // jet selection

    EffMinimCuts [iS] = 0; // mass, charge, vertex, cosmic
    EffAccCuts   [iS] = 0; // pt, eta
    EffIDCuts    [iS] = 0; // tight muon id
    EffIsoCuts   [iS] = 0; // trk Iso
    EffTrigCuts  [iS] = 0; // iso mu 24_eta2p1
    EffJetSel    [iS] = 0; // jet selection
    EffBDTvbfCut [iS] = 0; // VBF BDT cut
    EffNonJetSel [iS] = 0; // jet selection
    EffDiPt10GeV [iS] = 0; // pt(mumu) > = 10 GeV/c 

    dEffMinimCuts [iS] = 0; // mass, charge, vertex, cosmic
    dEffAccCuts   [iS] = 0; // pt, eta
    dEffIDCuts    [iS] = 0; // tight muon id
    dEffIsoCuts   [iS] = 0; // trk Iso
    dEffTrigCuts  [iS] = 0; // iso mu 24_eta2p1
    dEffJetSel    [iS] = 0; // jet selection
    dEffBDTvbfCut [iS] = 0; // VBF BDT cut
    dEffNonJetSel [iS] = 0; // jet selection
    dEffDiPt10GeV [iS] = 0; // pt(mumu) > = 10 GeV/c 
  }

  //////////////////////////////////////////////////////////////////////
  //create txt file with fit output
  ofstream myfile ("EffInfo.txt");
  std::string effFileName = outputFileName;
  effFileName.append(".txt");
  //ofstream myfileSubCat ("EffInfoSubCat.txt");
  ofstream myfileSubCat (effFileName.c_str());
  //////////////////////////////////////////////////////////////////////


  unsigned nEvents = tree->GetEntries();
  unsigned nSubEven = unsigned(min(int(maxEvents),int(nEvents))/Nbin);
  cout << "nEvents: " << nEvents  << " and size of nSubEvents = " << nSubEven << endl;

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

    //for efficiency:
    //find sub sample
    int IsampleGen = -1;
    if(Flag_subCat != 1){
      for(unsigned iS = 0; iS<Nbin ;iS++){
        if(i >= nSubEven*iS && i < nSubEven*(iS+1)){
          IsampleGen = iS;
        }
        if(iS == (Nbin -1) && i >= nSubEven*iS){
          IsampleGen = iS;
        }
      }
    }
    if(Flag_subCat == 1)IsampleGen = 0;
    //cout << "subcategory = " << IsampleGen << " for iEvent = " << i << " and Nbin = " << Nbin << " nSubEven = " << nSubEven << endl;  
    if (IsampleGen < 0) cout << " IsampleGen = -1 for Event = " << i << endl;
    if (IsampleGen < 0) continue;
    //end: find sub sample
    // calculate check for gen Mass
    float MassGenBoson = genHpostFSR.mass;
    //if (counterGenBoson>= 10000) continue;
    //if (counterGenBoson>= 50000) continue;
    // check for fhe doublets:
    pairOfDouble massPt(genHpostFSR.mass,genHpostFSR.pt);
    //if(MassGenBoson >= MinMassEff && MassGenBoson < MaxMassEff && uniqueGeneratedEvents.insert(massPt).second){
    if(uniqueGeneratedEvents.insert(massPt).second){ // no mass cuts at all for gen level only
      counterGenBoson[IsampleGen]++;
      //hMassGenBoson ->Fill(MassGenBoson);
    }

    if(reco1.pt<0. || reco2.pt<0.)
      continue;

    double weight = 1.0;
    double weightMuonEffUp   = 1.0;
    double weightMuonEffDown = 1.0;
    if (isSignal)
      {
        double randForSF = randomForSF.Rndm();
        if (runPeriod == "7TeV"){
          weight            *= weightFromSF_2011(randForSF,reco1,reco2, 0.   , 0.   ,0.    );
          weightMuonEffUp   *= weightFromSF_2011(randForSF,reco1,reco2, 0.005, 0.002, 0.002);
          weightMuonEffDown *= weightFromSF_2011(randForSF,reco1,reco2,-0.005,-0.002,-0.002);
        }else{
          weight            *= weightFromSF(randForSF,reco1,reco2, 0.   , 0.   , 0.   );
          weightMuonEffUp   *= weightFromSF(randForSF,reco1,reco2, 0.005, 0.002, 0.002);
          weightMuonEffDown *= weightFromSF(randForSF,reco1,reco2,-0.005,-0.002,-0.002);
        }
      }



    _dimuonMass_noMuscle=recoCandMass;
    _dimuonPt_noMuscle=recoCandPt;
    if (reco1.pt>reco2.pt)
      {
        _muonLead_pt_noMuscle=reco1.pt;
        _muonSub_pt_noMuscle=reco2.pt;
        _muonLead_passPFRelIso_noMuscle  = getPFRelIso (reco1) < 0.12 ? 1 : 0;
        _muonLead_passTrkRelIso_noMuscle = getTrkRelIso(reco1) < 0.10 ? 1 : 0;
        _muonSub_passPFRelIso_noMuscle  = getPFRelIso (reco2) < 0.12 ? 1 : 0;
        _muonSub_passTrkRelIso_noMuscle = getTrkRelIso(reco2) < 0.10 ? 1 : 0;
      }
    else
      {
        _muonLead_pt_noMuscle=reco2.pt;
        _muonSub_pt_noMuscle=reco1.pt;
        _muonLead_passPFRelIso_noMuscle  = getPFRelIso (reco2) < 0.12 ? 1 : 0;
        _muonLead_passTrkRelIso_noMuscle = getTrkRelIso(reco2) < 0.10 ? 1 : 0;
        _muonSub_passPFRelIso_noMuscle  = getPFRelIso (reco1) < 0.12 ? 1 : 0;
        _muonSub_passTrkRelIso_noMuscle = getTrkRelIso(reco1) < 0.10 ? 1 : 0;
      }

    float recoCandMassOrig = recoCandMass;
    TLorentzVector reco1Vec;
    TLorentzVector reco2Vec;
    reco1Vec.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,MASS_MUON);
    reco2Vec.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,MASS_MUON);

    //Matching gen muons to reco muons
    TLorentzVector reco1GenVec;
    TLorentzVector reco2GenVec;
    float recoGenDR1 = -999.;
    float recoGenDR2 = -999.;
    float recoOrigRes1 = -999.;
    float recoOrigRes2 = -999.;
    if (reco1GenPostFSR.pt > 0.0 && reco2GenPostFSR.pt > 0.0)
      {
        reco1GenVec.SetPtEtaPhiM(reco1GenPostFSR.pt,reco1GenPostFSR.eta,reco1GenPostFSR.phi,MASS_MUON);
        reco2GenVec.SetPtEtaPhiM(reco2GenPostFSR.pt,reco2GenPostFSR.eta,reco2GenPostFSR.phi,MASS_MUON);
        float deltaR_t1r1 = reco1GenVec.DeltaR(reco1Vec);
        float deltaR_t1r2 = reco1GenVec.DeltaR(reco2Vec);
        float deltaR_t2r1 = reco2GenVec.DeltaR(reco1Vec);
        float deltaR_t2r2 = reco2GenVec.DeltaR(reco2Vec);
        _TrackInfo tmpTrue1 = reco1GenPostFSR;
        _TrackInfo tmpTrue2 = reco2GenPostFSR;
        TLorentzVector tmpVec1 = reco1GenVec;
        TLorentzVector tmpVec2 = reco2GenVec;
        if (deltaR_t1r1 > deltaR_t1r2)
          {
            tmpTrue1 = reco2GenPostFSR;
            tmpVec1 = reco2GenVec;
          }
        if (deltaR_t2r2 > deltaR_t2r1)
          {
            tmpTrue2 = reco1GenPostFSR;
            tmpVec2 = reco1GenVec;
          }
        reco1GenPostFSR=tmpTrue1;
        reco2GenPostFSR=tmpTrue2;
        reco1GenVec = tmpVec1;
        reco2GenVec = tmpVec2;
        recoGenDR1 = reco1GenVec.DeltaR(reco1Vec);
        recoGenDR2 = reco2GenVec.DeltaR(reco2Vec);
        recoOrigRes1 = (reco1.pt-reco1GenPostFSR.pt)/reco1.pt;
        recoOrigRes2 = (reco2.pt-reco2GenPostFSR.pt)/reco2.pt;
      }
    // reject fake muon matches to gen level for signal only : very loose rejection DR < 0.5 and Resolution < 40%!
    //if(isSignal){
    //   if(recoGenDR1 > 0.5 || recoGenDR2 > 0.5) continue;
    //   if(fabs(recoOrigRes1) > 0.4 || fabs(recoOrigRes2) > 0.4) continue;
    //}
#ifdef MUSCLEFIT
    TLorentzVector reco1Cor;
    TLorentzVector reco2Cor;
    reco1Cor.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,MASS_MUON);
    reco2Cor.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,MASS_MUON);
    mfCorr->applyPtCorrection(reco1Cor,reco1.charge);
    mfCorr->applyPtCorrection(reco2Cor,reco2.charge);

    if (!isData && runPeriod=="8TeV"){
      mfCorr->applyPtSmearing(reco1Cor,reco1.charge);
      mfCorr->applyPtSmearing(reco2Cor,reco2.charge);
    }

    TLorentzVector diMuonCor = reco1Cor + reco2Cor;
    // recalculate ptErr which come from curverture of magnetic field:
    if(reco1.pt > 0)reco1.ptErr = reco1.ptErr*reco1Cor.Pt()/reco1.pt;
    if(reco2.pt > 0)reco2.ptErr = reco2.ptErr*reco2Cor.Pt()/reco2.pt;
    // fill corr  pT:
    reco1.pt = reco1Cor.Pt();
    reco2.pt = reco2Cor.Pt();
    recoCandMass = diMuonCor.M();
    recoCandPt = diMuonCor.Pt();
    recoCandY = diMuonCor.Rapidity();
    recoCandPhi = diMuonCor.Phi();
    reco1Vec = reco1Cor;
    reco2Vec = reco2Cor;
#endif
#ifdef ROCHESTER
    TLorentzVector reco1Cor;
    TLorentzVector reco2Cor;
    reco1Cor.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,MASS_MUON);
    reco2Cor.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,MASS_MUON);
    float rochesterError=1.0; //1.0 if you don't care
    if (runPeriod == "7TeV")
      {
        if (isData)
          {
            rCorr11->momcor_data(reco1Cor,reco1.charge,0,rochesterRun,rochesterError);
            rCorr11->momcor_data(reco2Cor,reco2.charge,0,rochesterRun,rochesterError);
          }
        else
          {
            rCorr11->momcor_mc(reco1Cor,reco1.charge,0,rochesterRun,rochesterError);
            rCorr11->momcor_mc(reco2Cor,reco2.charge,0,rochesterRun,rochesterError);
          }
      }
    else
      {
        if (isData)
          {
            rCorr12->momcor_data(reco1Cor,reco1.charge,0,rochesterRun,rochesterError);
            rCorr12->momcor_data(reco2Cor,reco2.charge,0,rochesterRun,rochesterError);
          }
        else
          {
            rCorr12->momcor_mc(reco1Cor,reco1.charge,0,rochesterRun,rochesterError);
            rCorr12->momcor_mc(reco2Cor,reco2.charge,0,rochesterRun,rochesterError);
          }
      }
    TLorentzVector diMuonCor = reco1Cor + reco2Cor;
    reco1.pt = reco1Cor.Pt();
    reco2.pt = reco2Cor.Pt();
    recoCandMass = diMuonCor.M();
    recoCandPt = diMuonCor.Pt();
    recoCandY = diMuonCor.Rapidity();
    recoCandPhi = diMuonCor.Phi();
    reco1Vec = recoCor1;
    reco2Vec = recoCor2;
#endif
#ifdef ROCHESTER2012ReReco
    TLorentzVector reco1Cor;
    TLorentzVector reco2Cor;
    reco1Cor.SetPtEtaPhiM(reco1.pt,reco1.eta,reco1.phi,MASS_MUON);
    reco2Cor.SetPtEtaPhiM(reco2.pt,reco2.eta,reco2.phi,MASS_MUON);
    float rochesterError=1.0; //1.0 if you don't care
    if (isData)
      {
        rCorr12jan22->momcor_data(reco1Cor,reco1.charge,0,rochesterError);
        rCorr12jan22->momcor_data(reco2Cor,reco2.charge,0,rochesterError);
      }
    else
      {
        rCorr12jan22->momcor_mc(reco1Cor,reco1.charge,0,rochesterError);
        rCorr12jan22->momcor_mc(reco2Cor,reco2.charge,0,rochesterError);
      }
    TLorentzVector diMuonCor = reco1Cor + reco2Cor;
    reco1.pt = reco1Cor.Pt();
    reco2.pt = reco2Cor.Pt();
    recoCandMass = diMuonCor.M();
    recoCandPt = diMuonCor.Pt();
    recoCandY = diMuonCor.Rapidity();
    recoCandPhi = diMuonCor.Phi();
    reco1Vec = reco1Cor;
    reco2Vec = reco2Cor;
#endif

    TLorentzVector recoCandVec = reco1Vec+reco2Vec;

    float mDiMuResSFUp   = recoCandMass;
    float mDiMuResSFDown = recoCandMass;

#ifdef SMEARING
    if(isSignal) // smear only signal because it has muons from higgs 
      {
        if(reco1GenPostFSR.pt<0.)
          cout << "Muon 1 Post FSR not valid!\n";
        if(reco2GenPostFSR.pt<0.)
          cout << "Muon 2 Post FSR not valid!\n";
        float ptReco1 = -1.;
        float ptReco2 = -1.;
        if(runPeriod == "7TeV")
          {
            ptReco1 = smearPT2011 -> PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR);
            ptReco2 = smearPT2011 -> PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR);
          }
        else
          {
            ptReco1 = smearPT -> PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR);
            ptReco2 = smearPT -> PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR);
          }
        TLorentzVector reco1Vec;
        TLorentzVector reco2Vec;
        reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
        reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
        TLorentzVector diMuonVec = reco1Vec + reco2Vec;

        reco1.pt = ptReco1;
        reco2.pt = ptReco2;
        recoCandMass = diMuonVec.M();
        recoCandPt = diMuonVec.Pt();
        recoCandY = diMuonVec.Rapidity();
        recoCandPhi = diMuonVec.Phi();
      
        //Systematics Time For Muon Momentum Resolution
        if(runPeriod == "7TeV")
          {
            ptReco1 = smearPT2011->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR, 0.13);
            ptReco2 = smearPT2011->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR, 0.13);
          }
        else
          {
            ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR, 0.02);
            ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR, 0.02);
          }
        reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
        reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
        diMuonVec = reco1Vec + reco2Vec;
        mDiMuResSFUp = diMuonVec.M();
      
        if(runPeriod == "7TeV")
          {
            ptReco1 = smearPT2011->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR, -0.13);
            ptReco2 = smearPT2011->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR, -0.13);
          }
        else
          {
            ptReco1 = smearPT->PTsmear(reco1GenPostFSR.pt, reco1GenPostFSR.eta, reco1GenPostFSR.charge, reco1.pt, ISMEAR, -0.02);
            ptReco2 = smearPT->PTsmear(reco2GenPostFSR.pt, reco2GenPostFSR.eta, reco2GenPostFSR.charge, reco2.pt, ISMEAR, -0.02);
          }
        reco1Vec.SetPtEtaPhiM(ptReco1,reco1.eta,reco1.phi,MASS_MUON);
        reco2Vec.SetPtEtaPhiM(ptReco2,reco2.eta,reco2.phi,MASS_MUON);
        diMuonVec = reco1Vec + reco2Vec;
        mDiMuResSFDown = diMuonVec.M();
      }
    
#endif

    // Anna's Extra Checks on Muons; Moved To After Corrections
    //if (counterGenBoson <= 40000) continue;
    // additional selection cuts
    //cout << "TEST iEVENT =  " << i  << " reco1.pt = " << reco1.pt << endl;
    if (reco1.charge == reco2.charge) continue;
    if(reco1.eta < -900. || reco2.eta < -900) continue;//rejection fake in reco level


    // EFFTEST:
#ifdef EFFTEST 

    // 1st EFFTEST

    if(reco1.charge != reco2.charge && reco1.pt > 20 && reco2.pt > 20 && fabs(reco1.eta)<2.4 && recoCandMass > 50){
      
      cout << "info1: " << eventInfo.run <<  " " << eventInfo.lumi << " " << eventInfo.event << " "
           << reco1.pt << " " << reco1.eta << " " << reco1.phi << " "
           << reco2.pt << " " << reco2.eta << " " << reco2.phi << " "
           << endl;

      if(isKinTight_2012_noIso_noPF(reco1) && isKinTight_2012_noIso_noPF(reco2)){

        cout << "info2 tight noPFIso, noPF: " << eventInfo.run <<  " " << eventInfo.lumi << " " << eventInfo.event << " "
             << reco1.pt << " " << reco1.eta << " " << reco1.phi << " "
             << reco2.pt << " " << reco2.eta << " " << reco2.phi << " "
             << endl;

      }
      if( (isKinTight_2012_noIso_noPF(reco1) && isKinTight_2012_noIso_noPF(reco2))
          && reco1.isPFMuon && reco1.isPFMuon  ){
        cout << "info3 tight noPFIso: " << eventInfo.run <<  " " << eventInfo.lumi << " " << eventInfo.event << " "
             << reco1.pt << " " << reco1.eta << " " << reco1.phi << " "
             << reco2.pt << " " << reco2.eta << " " << reco2.phi << " "
             << endl;

      }
    }
    // END EFFTEST
#endif

    fillMuonHist(hists.countsHist2, reco1, reco2);
    //printStationMiss(reco1,reco2,eventInfo,testString,testCounter);

    mva.resetValues();
    mva.mDiMu = recoCandMass;
    mva.RelMassRes = 0.;
    mva.RelMassResCov = 0.;
    //if(recoCandMassResCov < 0.0001) cout << "recoCandMassResCov = " << recoCandMassResCov << endl;
    if(recoCandMass > 0) 
      mva.RelMassRes = GetMassRes(reco1,reco2)/recoCandMass;
    if(recoCandMassOrig > 0) 
      mva.RelMassResCov = recoCandMassResCov/recoCandMassOrig;


    bool inBlindWindow = mva.mDiMu < maxBlind && mva.mDiMu > minBlind;

    bool blind = false;
  
#ifdef BLIND
    if (inBlindWindow && isData)
      blind = true;
#endif

#ifdef PUREWEIGHT
    if (!isData)
      {
        weight            *= lumiWeights.weight(nPU);
        weightMuonEffUp   *= lumiWeights.weight(nPU);
        weightMuonEffDown *= lumiWeights.weight(nPU);

      }
#endif

    hists.countsHist->Fill(0.0, weight);

    // make order in pt for muons
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
    // define geometrical categories: 
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

    // for efficiency:
    //find sub sample
    int Isample = -1;
    if(Flag_subCat != 1){
      for(unsigned iS = 0; iS<Nbin ;iS++){
        if(i >= nSubEven*iS && i < nSubEven*(iS+1)){
          Isample = iS;
        }
        if(iS == (Nbin -1) && i >= nSubEven*iS){
          Isample = iS;
        }
      }
    }
    if(Flag_subCat == 1){
      if(isBB) Isample = 0;
      if(isBO) Isample = 1;
      if(isBE) Isample = 2;
      if(isOO) Isample = 3;
      if(isOE) Isample = 4;
      if(isEE) Isample = 5;
    }
    //cout << "subcategory = " << Isample << " for iEvent = " << i << " and Nbin = " << Nbin << " nSubEven = " << nSubEven << endl;  
    if (Isample < 0) cout << " Isample = -1 for Event = " << i << endl;
    if (Isample < 0) continue;
    //end: find sub sample

    bool IFMinimCuts = mva.mDiMu >= MinMassEff && mva.mDiMu <= MaxMassEff;
    bool IFAccCuts   = reco1.pt >= 25 && reco2.pt >= 25 && fabs(reco1.eta) <= 2.1 && fabs(reco2.eta) <= 2.1; 
    bool IFIDCuts    = ((*muonIdFuncPtr)(reco1)) && ((*muonIdFuncPtr)(reco2)); // NOT contain ISO anymore
    bool IFPFIsoCuts   = getPFRelIso(reco1) <= 0.12 && getPFRelIso(reco2) <= 0.12;
    bool IFTrigCuts  = isHltMatched(reco1,reco2,allowedHLTPaths);
    bool IFAllTrigCuts  = IFMinimCuts && IFAccCuts && IFIDCuts && IFPFIsoCuts && IFTrigCuts;

    if(IFMinimCuts)                                          counterMinimCuts[Isample]++;
    if(IFMinimCuts && IFAccCuts)                             counterAccCuts[Isample]++;
    if(IFMinimCuts && IFAccCuts && IFIDCuts)                 counterIDCuts[Isample]++;
    if(IFMinimCuts && IFAccCuts && IFIDCuts && IFPFIsoCuts)  counterIsoCuts[Isample]++;
    if(IFAllTrigCuts)                                        counterTrigCuts[Isample]++; 

    // this selection does NOT contain ISO anymore!!!
    if (!((*muonIdFuncPtr)(reco1)) || !((*muonIdFuncPtr)(reco2)))
      continue;

    hists.countsHist->Fill(1.0, weight);

    // moved below
    //if (!isHltMatched(reco1,reco2,allowedHLTPaths))
    //    continue;

    hists.countsHist->Fill(2.0, weight);

    if (reco1.charge*reco2.charge != -1)
      continue;

    hists.countsHist->Fill(3.0, weight);

    if (mva.mDiMu < minMmm || mva.mDiMu > maxMmm)
      continue;

    hists.countsHist->Fill(4.0, weight);


    mva.weight = weight;
    mva.met = met.pt;
    mva.nPU = nPU;

    mva.ptMu1=muon1.pt;
    mva.ptMu2=muon2.pt;
    mva.etaMu1=muon1.eta;
    mva.etaMu2=muon2.eta;
    mva.deltaEtaMuons=fabs(muon1.eta-muon2.eta);
    
    mva.relIsoMu1    = getPFRelIso (muon1);
    mva.relIsoMu2    = getPFRelIso (muon2);
    mva.trkRelIsoMu1 = getTrkRelIso(muon1);
    mva.trkRelIsoMu2 = getTrkRelIso(muon2);

    mva.ptDiMu = recoCandPt;
    mva.yDiMu = recoCandY;

    float mDiMuCalibUp = mva.mDiMu+calibSysSmear;
    float mDiMuCalibDown = mva.mDiMu-calibSysSmear;


    float mDiMuResUp   = smearMC(trueMass,recoCandMass,calib,resSmear+resSysSmear,random);
    float mDiMuResDown = smearMC(trueMass,recoCandMass,calib,resSmear-resSysSmear,random);

    bool inTrainingWindow = (mva.mDiMu < 160. && mva.mDiMu > 70.);
    
    //////////////////////////////////////////
    // MEKD

    _kd = -1.; 
    _sigME = -1.; 
    _bakME = -1.;
    _kdPdf = -1.; 
    _sigMEPdf = -1.; 
    _bakMEPdf = -1.;
#ifdef MEKD_STANDALONE
    int kdStatus = 0;
    kdStatus = mekd.getKD(reco1Vec, reco2Vec, 
                          reco1.charge, 
                          _kd, _sigME, _bakME);
    if (kdStatus != 0)
      cout << "Error: MEKD: " << _kd << " Status: "<<kdStatus 
           << " sigME: " << _sigME << " bakME: "<< _bakME<<endl;

    kdStatus = mekdPdf.getKD(reco1Vec, reco2Vec, 
                             reco1.charge, 
                             _kdPdf, _sigMEPdf, _bakMEPdf);
    if (kdStatus != 0)
      cout << "Error: MEKD w/ PDF: " << _kdPdf << " Status: "<<kdStatus 
           << " sigME: " << _sigMEPdf << " bakME: "<< _bakMEPdf<<endl;
#endif

    //////////////////////////////////////////
    //Computing CosTheta*

    TLorentzVector pMuon1;
    TLorentzVector pMuon2;
    pMuon1.SetPtEtaPhiM(muon1.pt,muon1.eta,muon1.phi,MASS_MUON);
    pMuon2.SetPtEtaPhiM(muon2.pt,muon2.eta,muon2.phi,MASS_MUON);
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

    if (!blind)
      {
        hists.mDiMu->Fill(mva.mDiMu, weight);
        hists.yVmDiMu->Fill(mva.mDiMu,fabs(mva.yDiMu), weight);
        hists.ptVmDiMu->Fill(mva.mDiMu,mva.ptDiMu, weight);
        hists.phiVmDiMu->Fill(mva.mDiMu,recoCandPhi, weight);

        hists.RelMassRes->Fill(mva.RelMassRes, weight);
        hists.RelMassResCov->Fill(mva.RelMassResCov, weight);
      }

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
    hists.trkRelIsoMu1->Fill(mva.trkRelIsoMu1, weight);
    hists.trkRelIsoMu2->Fill(mva.trkRelIsoMu2, weight);

    hists.nPU->Fill(nPU, weight);
    hists.nVtx->Fill(mva.nVtx, weight);
    hists.met->Fill(met.pt, weight);
    hists.weight->Fill(weight);

    // Jet Part
    _PFJetInfo originalJets = jets;
    float jetPtCut = 30.;
    float jetPtCutC = 30.;
    int jetPUIDCut = 4; // >=    tight = 7, medium = 6, loose = 4
    std::vector<unsigned> goodJetIndices;
    for(unsigned iJet=0; (iJet < jets.nJets && iJet < 10);iJet++)
      {
        if (jets.genPt[iJet]>0.0 && jets.pt[iJet]>15.)
          jets.pt[iJet] = jerCorr(jets.pt[iJet],jets.genPt[iJet],jets.eta[iJet]);
        bool goodPt = jets.pt[iJet]>jetPtCut || (fabs(jets.eta[iJet < 2.4]) && jets.pt[iJet]>jetPtCutC);
        bool goodPUID = puJetFullId[iJet] >= jetPUIDCut;
        if (goodPt && goodPUID)
          {
            mva.nJets++;
            mva.ht += jets.pt[iJet];
            goodJetIndices.push_back(iJet);
          }
      }

    if (mva.nJets>=1)
    {
      unsigned iJet1 = goodJetIndices[0];
      TLorentzVector pJet1;
      pJet1.SetXYZM(jets.px[iJet1],jets.py[iJet1],jets.pz[iJet1],jets.mass[iJet1]);

      mva.ptJet1 = pJet1.Pt();
      mva.etaJet1 = pJet1.Eta();
      hists.ptJet1->Fill(mva.ptJet1, weight);
      hists.etaJet1->Fill(mva.etaJet1, weight);
      _jetLead_flav = jets.partonFlavour[iJet1];

      mva.puJetIDSimpleDiscJet1 = puJetSimpleDisc[iJet1];
      mva.puJetIDFullDiscJet1 = puJetFullDisc[iJet1];
      mva.puJetIDSimpleJet1 = (int) puJetSimpleId[iJet1];
      mva.puJetIDFullJet1 = (int) puJetFullId[iJet1];
      hists.puJetIDSimpleDiscJet1->Fill(mva.puJetIDSimpleDiscJet1,weight);
      hists.puJetIDSimpleJet1->Fill(mva.puJetIDSimpleJet1,weight);

      mva.deltaPhiHJ1 = pJet1.DeltaPhi(diMuon);
      hists.deltaPhiHJ1->Fill(mva.deltaPhiHJ1, weight);
      if (mva.nJets < 2)
      {
        unsigned iJet1 = goodJetIndices[0];
        TLorentzVector pJet1;
        pJet1.SetXYZM(jets.px[iJet1],jets.py[iJet1],jets.pz[iJet1],jets.mass[iJet1]);

        mva.ptJet1 = pJet1.Pt();
        mva.etaJet1 = pJet1.Eta();
        hists.ptJet1->Fill(mva.ptJet1, weight);
        hists.etaJet1->Fill(mva.etaJet1, weight);

        mva.puJetIDSimpleDiscJet1 = puJetSimpleDisc[iJet1];
        mva.puJetIDFullDiscJet1 = puJetFullDisc[iJet1];
        mva.puJetIDSimpleJet1 = (int) puJetSimpleId[iJet1];
        mva.puJetIDFullJet1 = (int) puJetFullId[iJet1];
        hists.puJetIDSimpleDiscJet1->Fill(mva.puJetIDSimpleDiscJet1,weight);
        hists.puJetIDSimpleJet1->Fill(mva.puJetIDSimpleJet1,weight);

        mva.deltaPhiHJ1 = pJet1.DeltaPhi(diMuon);
        hists.deltaPhiHJ1->Fill(mva.deltaPhiHJ1, weight);
        if (mva.nJets < 2)
          {
            mva.ptmiss = (pJet1+recoCandVec).Pt();
          }
      }
    }

    if(mva.nJets>=2)
    {
      unsigned iJet1 = goodJetIndices[0];
      unsigned iJet2 = goodJetIndices[1];
      TLorentzVector pJet1;
      TLorentzVector pJet2;
      pJet1.SetXYZM(jets.px[iJet1],jets.py[iJet1],jets.pz[iJet1],jets.mass[iJet1]);
      pJet2.SetXYZM(jets.px[iJet2],jets.py[iJet2],jets.pz[iJet2],jets.mass[iJet2]);
      TLorentzVector diJet = pJet1+pJet2;

      double dEtaJets = fabs(jets.eta[iJet1]-jets.eta[iJet2]);
      double etaJetProduct = jets.eta[iJet1]*jets.eta[iJet2];
      mva.deltaPhiJets = pJet1.DeltaPhi(pJet2);
      mva.deltaRJets = pJet1.DeltaR(pJet2);

      // Seeing if there are jets in the rapidity gap
      float etaMax = jets.eta[iJet1];
      float etaMin = 9999999.0;
      if(etaMax < jets.eta[iJet2])
      {
          etaMax = jets.eta[iJet2];
          etaMin = jets.eta[iJet1];
      }
      else
      {
          etaMin = jets.eta[iJet2];
      }
      bool jetInRapidityGap=false;
      for(std::vector<unsigned>::const_iterator iGoodJet=goodJetIndices.begin();
                                      iGoodJet != goodJetIndices.end();iGoodJet++)
      {
          if(jets.eta[*iGoodJet] < etaMax && jets.eta[*iGoodJet] > etaMin)
          {
            jetInRapidityGap = true;
            mva.nJetsInRapidityGap++;
            mva.htInRapidityGap += jets.pt[*iGoodJet];
          }
      }

      _jetSub_flav = jets.partonFlavour[iJet2];

      mva.puJetIDSimpleDiscJet2 = puJetSimpleDisc[iJet2];
      mva.puJetIDFullDiscJet2 = puJetFullDisc[iJet2];
      mva.puJetIDSimpleJet2 = (int) puJetSimpleId[iJet2];
      mva.puJetIDFullJet2 = (int) puJetFullId[iJet2];

      if (mva.nJets>=3)
      {
        unsigned iJet3 = goodJetIndices[2];
        mva.puJetIDSimpleDiscJet3 = puJetSimpleDisc[iJet3];
        mva.puJetIDFullDiscJet3 = puJetFullDisc[iJet3];
        mva.puJetIDSimpleJet3 = (int) puJetSimpleId[iJet3];
        mva.puJetIDFullJet3 = (int) puJetFullId[iJet3];
        hists.puJetIDSimpleDiscJet3->Fill(mva.puJetIDSimpleDiscJet3,weight);
        hists.puJetIDSimpleJet3->Fill(mva.puJetIDSimpleJet3,weight);
      }

      mva.mDiJet = diJet.M();
      mva.yDiJet = diJet.Rapidity();
      mva.ptDiJet = diJet.Pt();
      mva.ptJet2 = pJet2.Pt();
      mva.etaJet2 = pJet2.Eta();
      mva.productEtaJets = etaJetProduct;
      mva.deltaEtaJets = dEtaJets;
      mva.ptmiss = (diJet+recoCandVec).Pt();

      hists.mDiJet->Fill(mva.mDiJet, weight);
      hists.ptDiJet->Fill(mva.ptDiJet, weight);
      hists.yDiJet->Fill(mva.yDiJet, weight);
      hists.ptJet2->Fill(mva.ptJet2, weight);
      hists.etaJet2->Fill(mva.etaJet2, weight);
      hists.deltaEtaJets->Fill(mva.deltaEtaJets, weight);
      hists.deltaPhiJets->Fill(mva.deltaPhiJets, weight);
      hists.deltaRJets->Fill(mva.deltaRJets, weight);
      hists.nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, weight);
      hists.htInRapidityGap->Fill(mva.htInRapidityGap, weight);
      hists.nJets->Fill(mva.nJets, weight);
      hists.ht->Fill(mva.ht, weight);
      hists.ptmiss->Fill(mva.ptmiss, weight);

      hists.puJetIDSimpleDiscJet2->Fill(mva.puJetIDSimpleDiscJet2,weight);

      hists.puJetIDSimpleJet2->Fill(mva.puJetIDSimpleJet2,weight);
    }

    /////////////////////////////////////////////
    // JES Uncertainties
  
    _nJets_JESUp = 0.0;
    _ptMiss_JESUp = -1.0;
    _deltaEtaJets_JESUp = -1.0;
    _dijetMass_JESUp = -1.0;
  
    _jetLead_pt_JESUp = -1.0;          
    _jetLead_eta_JESUp = -999.0;         
    _jetSub_pt_JESUp = -1.0;          
    _jetSub_eta_JESUp = -999.0;         
  
    _nJets_JESDown = 0.0;
    _ptMiss_JESDown = -1.0;
    _deltaEtaJets_JESDown = -1.0;
    _dijetMass_JESDown = -1.0;
  
    _jetLead_pt_JESDown = -1.0;          
    _jetLead_eta_JESDown = -999.0;         
    _jetSub_pt_JESDown = -1.0;          
    _jetSub_eta_JESDown = -999.0;         
  
    // JES Up
    _PFJetInfo jetsJESUp = jets;
    goodJetIndices.clear();
    for(unsigned iJet=0; (iJet < jetsJESUp.nJets && iJet < 10);iJet++)
      {
        if (jetsJESUp.jecUnc[iJet]>0.0 && jetsJESUp.jecUnc[iJet]<1.)
          {
            jetsJESUp.pt[iJet] += jetsJESUp.jecUnc[iJet]*jetsJESUp.pt[iJet];
          }
        bool goodPt = jetsJESUp.pt[iJet]>jetPtCut || (fabs(jetsJESUp.eta[iJet < 2.4]) && jetsJESUp.pt[iJet]>jetPtCutC);
        bool goodPUID = puJetFullId[iJet] >= jetPUIDCut;
        if (goodPt && goodPUID)
          {
            _nJets_JESUp++;
            goodJetIndices.push_back(iJet);
          }
      }

    if (goodJetIndices.size()>=1)
      {
        unsigned iJet1 = goodJetIndices[0];
        _jetLead_pt_JESUp = jetsJESUp.pt[iJet1];
        _jetLead_eta_JESUp = jetsJESUp.eta[iJet1];
      }

    if(goodJetIndices.size()>=2)
      {
        unsigned iJet1 = goodJetIndices[0];
        unsigned iJet2 = goodJetIndices[1];
        TLorentzVector pJet1;
        TLorentzVector pJet2;
        pJet1.SetXYZM(jetsJESUp.px[iJet1],jetsJESUp.py[iJet1],jetsJESUp.pz[iJet1],jetsJESUp.mass[iJet1]);
        pJet2.SetXYZM(jetsJESUp.px[iJet2],jetsJESUp.py[iJet2],jetsJESUp.pz[iJet2],jetsJESUp.mass[iJet2]);
        TLorentzVector diJet = pJet1+pJet2;

        double dEtaJets = fabs(jetsJESUp.eta[iJet1]-jetsJESUp.eta[iJet2]);
        double etaJetProduct = jetsJESUp.eta[iJet1]*jetsJESUp.eta[iJet2];

        _dijetMass_JESUp = diJet.M();
        _jetSub_pt_JESUp = pJet2.Pt();
        _jetSub_eta_JESUp = pJet2.Eta();
        _deltaEtaJets_JESUp = dEtaJets;
        _ptMiss_JESUp = (diJet+recoCandVec).Pt();
      }

    // JES Down
    _PFJetInfo jetsJESDown = jets;
    goodJetIndices.clear();
    for(unsigned iJet=0; (iJet < jetsJESDown.nJets && iJet < 10);iJet++)
      {
        if (jetsJESDown.jecUnc[iJet]>0.0 && jetsJESDown.jecUnc[iJet]<1.)
          {
            jetsJESDown.pt[iJet] -= jetsJESDown.jecUnc[iJet]*jetsJESDown.pt[iJet];
          }
        bool goodPt = jetsJESDown.pt[iJet]>jetPtCut || (fabs(jetsJESDown.eta[iJet < 2.4]) && jetsJESDown.pt[iJet]>jetPtCutC);
        bool goodPUID = puJetFullId[iJet] >= jetPUIDCut;
        if (goodPt && goodPUID)
          {
            _nJets_JESDown++;
            goodJetIndices.push_back(iJet);
          }
      }

    if (goodJetIndices.size()>=1)
      {
        unsigned iJet1 = goodJetIndices[0];
        _jetLead_pt_JESDown = jetsJESDown.pt[iJet1];
        _jetLead_eta_JESDown = jetsJESDown.eta[iJet1];
      }

    if(goodJetIndices.size()>=2)
      {
        unsigned iJet1 = goodJetIndices[0];
        unsigned iJet2 = goodJetIndices[1];
        TLorentzVector pJet1;
        TLorentzVector pJet2;
        pJet1.SetXYZM(jetsJESDown.px[iJet1],jetsJESDown.py[iJet1],jetsJESDown.pz[iJet1],jetsJESDown.mass[iJet1]);
        pJet2.SetXYZM(jetsJESDown.px[iJet2],jetsJESDown.py[iJet2],jetsJESDown.pz[iJet2],jetsJESDown.mass[iJet2]);
        TLorentzVector diJet = pJet1+pJet2;

        double dEtaJets = fabs(jetsJESDown.eta[iJet1]-jetsJESDown.eta[iJet2]);
        double etaJetProduct = jetsJESDown.eta[iJet1]*jetsJESDown.eta[iJet2];

        _dijetMass_JESDown = diJet.M();
        _jetSub_pt_JESDown = pJet2.Pt();
        _jetSub_eta_JESDown = pJet2.Eta();
        _deltaEtaJets_JESDown = dEtaJets;
        _ptMiss_JESDown = (diJet+recoCandVec).Pt();
      }

    /////////////////////////////////////////////
    // JER Uncertainties
  
    _nJets_JERUp = 0.0;
    _ptMiss_JERUp = -1.0;
    _deltaEtaJets_JERUp = -1.0;
    _dijetMass_JERUp = -1.0;
  
    _jetLead_pt_JERUp = -1.0;          
    _jetLead_eta_JERUp = -999.0;         
    _jetSub_pt_JERUp = -1.0;          
    _jetSub_eta_JERUp = -999.0;         
  
    _nJets_JERDown = 0.0;
    _ptMiss_JERDown = -1.0;
    _deltaEtaJets_JERDown = -1.0;
    _dijetMass_JERDown = -1.0;
  
    _jetLead_pt_JERDown = -1.0;          
    _jetLead_eta_JERDown = -999.0;         
    _jetSub_pt_JERDown = -1.0;          
    _jetSub_eta_JERDown = -999.0;         
  
    // JER Up
    _PFJetInfo jetsJERUp = originalJets;
    goodJetIndices.clear();
    for(unsigned iJet=0; (iJet < jetsJERUp.nJets && iJet < 10);iJet++)
      {
        if (jetsJERUp.genPt[iJet]>0.0 && jetsJERUp.pt[iJet]>15.)
          {
            jetsJERUp.pt[iJet] = corrPtUp(jetsJERUp.pt[iJet],jetsJERUp.genPt[iJet],jetsJERUp.eta[iJet]);
          }
        bool goodPt = jetsJERUp.pt[iJet]>jetPtCut || (fabs(jetsJERUp.eta[iJet < 2.4]) && jetsJERUp.pt[iJet]>jetPtCutC);
        bool goodPUID = puJetFullId[iJet] >= jetPUIDCut;
        if (goodPt && goodPUID)
          {
            _nJets_JERUp++;
            goodJetIndices.push_back(iJet);
          }
      }

    if (goodJetIndices.size()>=1)
      {
        unsigned iJet1 = goodJetIndices[0];
        _jetLead_pt_JERUp = jetsJERUp.pt[iJet1];
        _jetLead_eta_JERUp = jetsJERUp.eta[iJet1];
      }

    if(goodJetIndices.size()>=2)
      {
        unsigned iJet1 = goodJetIndices[0];
        unsigned iJet2 = goodJetIndices[1];
        TLorentzVector pJet1;
        TLorentzVector pJet2;
        pJet1.SetXYZM(jetsJERUp.px[iJet1],jetsJERUp.py[iJet1],jetsJERUp.pz[iJet1],jetsJERUp.mass[iJet1]);
        pJet2.SetXYZM(jetsJERUp.px[iJet2],jetsJERUp.py[iJet2],jetsJERUp.pz[iJet2],jetsJERUp.mass[iJet2]);
        TLorentzVector diJet = pJet1+pJet2;

        double dEtaJets = fabs(jetsJERUp.eta[iJet1]-jetsJERUp.eta[iJet2]);
        double etaJetProduct = jetsJERUp.eta[iJet1]*jetsJERUp.eta[iJet2];

        _dijetMass_JERUp = diJet.M();
        _jetSub_pt_JERUp = pJet2.Pt();
        _jetSub_eta_JERUp = pJet2.Eta();
        _deltaEtaJets_JERUp = dEtaJets;
        _ptMiss_JERUp = (diJet+recoCandVec).Pt();
      }

    // JER Down
    _PFJetInfo jetsJERDown = originalJets;
    goodJetIndices.clear();
    for(unsigned iJet=0; (iJet < jetsJERDown.nJets && iJet < 10);iJet++)
      {
        if (jetsJERDown.genPt[iJet]>0.0 && jetsJERDown.pt[iJet]>15.)
          {
            jetsJERDown.pt[iJet] = corrPtDown(jetsJERDown.pt[iJet],jetsJERDown.genPt[iJet],jetsJERDown.eta[iJet]);
          }
        bool goodPt = jetsJERDown.pt[iJet]>jetPtCut || (fabs(jetsJERDown.eta[iJet < 2.4]) && jetsJERDown.pt[iJet]>jetPtCutC);
        bool goodPUID = puJetFullId[iJet] >= jetPUIDCut;
        if (goodPt && goodPUID)
          {
            _nJets_JERDown++;
            goodJetIndices.push_back(iJet);
          }
      }

    if (goodJetIndices.size()>=1)
      {
        unsigned iJet1 = goodJetIndices[0];
        _jetLead_pt_JERDown = jetsJERDown.pt[iJet1];
        _jetLead_eta_JERDown = jetsJERDown.eta[iJet1];
      }

    if(goodJetIndices.size()>=2)
      {
        unsigned iJet1 = goodJetIndices[0];
        unsigned iJet2 = goodJetIndices[1];
        TLorentzVector pJet1;
        TLorentzVector pJet2;
        pJet1.SetXYZM(jetsJERDown.px[iJet1],jetsJERDown.py[iJet1],jetsJERDown.pz[iJet1],jetsJERDown.mass[iJet1]);
        pJet2.SetXYZM(jetsJERDown.px[iJet2],jetsJERDown.py[iJet2],jetsJERDown.pz[iJet2],jetsJERDown.mass[iJet2]);
        TLorentzVector diJet = pJet1+pJet2;

        double dEtaJets = fabs(jetsJERDown.eta[iJet1]-jetsJERDown.eta[iJet2]);
        double etaJetProduct = jetsJERDown.eta[iJet1]*jetsJERDown.eta[iJet2];

        _dijetMass_JERDown = diJet.M();
        _jetSub_pt_JERDown = pJet2.Pt();
        _jetSub_eta_JERDown = pJet2.Eta();
        _deltaEtaJets_JERDown = dEtaJets;
        _ptMiss_JERDown = (diJet+recoCandVec).Pt();
      }
  
  
    //HIG-12-007 PAS H->tautau
    //The VBF category requires at least two jets with pT > 30 GeV/c, |η1 − η2 | > 4.0 and
    //η1 · η2 < 0 (with η1 the pseudorapidty of the leading jet and η2 the pseudorapidity
    //of the subleading jet), and a di-jet invariant mass m12 > 400 GeV/c2 , with no other
    //jet with pT > 30 GeV/c in the rapidity region between the two jets.

    if (
        !(inBlindWindow && isData) && inTrainingWindow
        && isHltMatched(reco1,reco2,allowedHLTPaths)
        && mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12
        )
      mva.writeEvent();

    if (trainingTreeRun) //Skip Filling of histos when training Tree
      continue;

    bool vbfPreselection = mva.mDiJet>300.0 && mva.deltaEtaJets>3.0 && mva.productEtaJets<0.0;
    //if(vbfPreselection)
    //  std::cout << "VBF Preselected!!";
    mva.vbfPreselection = vbfPreselection;

    if(vbfPreselection)
      hists.countsHist->Fill(6);
    else
      hists.countsHist->Fill(5);
    
    // for efficiency:
    if(Flag_subCat != 1){
      if((!vbfPreselection) && IFAllTrigCuts ) counterNonJetSel[Isample]++;
      if ((!vbfPreselection) && mva.ptDiMu >= 10. && IFAllTrigCuts) counterDiPt10GeV[Isample]++;
      if(vbfPreselection && IFAllTrigCuts) counterJetSel[Isample]++;
    }
    if(Flag_subCat == 1){
      if((!vbfPreselection) && IFAllTrigCuts){
        counterNonJetSel[Isample]++;
        counterNonJetSel[6]++;
      }
      if ((!vbfPreselection) && mva.ptDiMu >= 10. && IFAllTrigCuts) {
        counterDiPt10GeV[Isample]++;
        counterDiPt10GeV[7]++;
        if(mva.RelMassRes < 0.009 && isBB) counterDiPt10GeV[8]++; //BB res1
        if(mva.RelMassRes >= 0.009 && isBB) counterDiPt10GeV[9]++; //BB res2 
        if(mva.RelMassRes < 0.011 && isBO) counterDiPt10GeV[10]++; //BO res1 
        if(mva.RelMassRes >= 0.011 && isBO) counterDiPt10GeV[11]++; //BO res2
      }
      if(vbfPreselection && IFAllTrigCuts){
        counterJetSel[Isample]++;
        counterJetSel[8]++;
        if(mva.mDiJet> 550.)counterJetSel[9]++;
        if(mva.mDiJet> 550. && mva.deltaEtaJets>3.5 && mva.ptmiss < 100)counterJetSel[10]++;
        if(mva.mDiJet> 500. && mva.deltaEtaJets>3.4 && mva.ptmiss < 25)counterJetSel[11]++;
      }
    }
    // 

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
    bool passIncBDTCut = mva.getMVAPassBDTCut(cfgNameInc);
    bool passVBFBDTCut = mva.getMVAPassBDTCut(cfgNameVBF);
    // for efficiency:
    if(IFAllTrigCuts && vbfPreselection && passVBFBDTCut && (Flag_subCat != 1)) counterBDTvbfCut[Isample]++;
    if(IFAllTrigCuts && vbfPreselection && passVBFBDTCut && (Flag_subCat == 1)){
      counterBDTvbfCut[Isample]++;
      counterBDTvbfCut[9]++;
    }
    //

    mva.bdtValInc = bdtValInc;
    mva.bdtValVBF = bdtValVBF;

    time_t timeStartFilling = time(NULL);

    if (!blind)
      {
        if(!vbfPreselection)
          {
            hists.BDTHistMuonOnly->Fill(bdtValInc, weight);
            hists.BDTHistMuonOnlyVMass->Fill(mva.mDiMu, bdtValInc, weight);
  
          }
        else
          {
            hists.BDTHistVBFVMass->Fill(mva.mDiMu, bdtValVBF, weight);
            hists.BDTHistVBF->Fill(bdtValVBF, weight);
          }
      }


    ///////////////////////////////////////////
    // Fill the outTree
    _eventInfo_run   = eventInfo.run  ;
    _eventInfo_lumi  = eventInfo.lumi ;
    _eventInfo_event = eventInfo.event;

    _nPU = mva.nPU;
    _puWeight = mva.weight;
    _puWeightMuonEffUp   = weightMuonEffUp;
    _puWeightMuonEffDown = weightMuonEffDown;
    _eventType = whichSelection(muon1,muon2,
                                runPeriod,
                                jets,
                                passIncBDTCut,
                                passVBFBDTCut);
    
    
    _dimuonMass   = mva.mDiMu;
    _dimuonMassResSFUp   = mDiMuResSFUp;
    _dimuonMassResSFDown = mDiMuResSFDown;
    _dimuonPt     = mva.ptDiMu;
    _dimuonY      = mva.yDiMu;
    _cosThetaStar = mva.cosThetaStar;

    _muonLead_charge = muon1.charge;
    _muonLead_pt     = muon1.pt;
    _muonLead_ptErr  = muon1.ptErr;
    _muonLead_eta    = muon1.eta;
    _muonLead_phi    = muon1.phi;
    _muonLead_passPFRelIso  = getPFRelIso (muon1) < 0.12 ? 1 : 0;
    _muonLead_passTrkRelIso = getTrkRelIso(muon1) < 0.10 ? 1 : 0;
    _muonLead_isHltMatched  = muon1.isHltMatched[0];

    _muonSub_charge = muon2.charge;
    _muonSub_pt     = muon2.pt;
    _muonSub_ptErr  = muon2.ptErr;
    _muonSub_eta    = muon2.eta;
    _muonSub_phi    = muon2.phi;
    _muonSub_passPFRelIso  = getPFRelIso (muon2) < 0.12 ? 1 : 0;
    _muonSub_passTrkRelIso = getTrkRelIso(muon2) < 0.10 ? 1 : 0;
    _muonSub_isHltMatched  = muon2.isHltMatched[0];
    
    _nJets = mva.nJets;
    _ptMiss = mva.ptmiss;
    _deltaEtaJets = mva.deltaEtaJets;
    _dijetMass       = mva.mDiJet;
    _dijetPt       = mva.ptDiJet;
    _dijetY       = mva.yDiJet;

    _jetLead_pt = mva.ptJet1;
    _jetLead_eta = mva.etaJet1;
    _jetLead_PUIDDisc = mva.puJetIDSimpleDiscJet1;
    _jetLead_PUIDFlag = mva.puJetIDSimpleJet1;
    _jetLead_FullPUIDDisc = mva.puJetIDFullDiscJet1;
    _jetLead_FullPUIDFlag = mva.puJetIDFullJet1;
    _jetLead_CutPUIDDisc = mva.puJetIDCutDiscJet1;

    _jetSub_pt = mva.ptJet2;
    _jetSub_eta = mva.etaJet2;
    _jetSub_PUIDDisc = mva.puJetIDSimpleDiscJet2;
    _jetSub_PUIDFlag = mva.puJetIDSimpleJet2;
    _jetSub_FullPUIDDisc = mva.puJetIDFullDiscJet2;
    _jetSub_FullPUIDFlag = mva.puJetIDFullJet2;
    _jetSub_CutPUIDDisc = mva.puJetIDCutDiscJet2;
    _bdtVBF = mva.bdtValVBF;
    _bdtNonVBF = mva.bdtValInc;

    // puid uncertainty weight.
    _puidUncWeight = 1.;
    if (fabs(mva.etaJet1)>2.4)
    {
      if (mva.ptJet1>40.)
        _puidUncWeight *= 1.02;
      else if(mva.ptJet1>30.)
        _puidUncWeight *= 1.05;
    }
    if (fabs(mva.etaJet2)>2.4)
    {
      if (mva.ptJet2>40.)
        _puidUncWeight *= 1.02;
      else if(mva.ptJet2>30.)
        _puidUncWeight *= 1.05;
    }

    if (_eventType != 0) _outTree -> Fill();

    ///////////////////////////////////////////

    // baseline++
    bool atLeastOneMuonFired = false;
    if ( reco1.pt>25 && isHltMatched(reco1,allowedHLTPaths) ) atLeastOneMuonFired = true;
    if ( reco2.pt>25 && isHltMatched(reco2,allowedHLTPaths) ) atLeastOneMuonFired = true;

    if ( atLeastOneMuonFired ) {

      //std::cout << "test" << std::endl;
      bool CutMassISO = IFMinimCuts && mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12; // mass cut + PF iso
      
      bool Jet2PtCuts = false;

      if ( mva.ptJet1 > 40. &&
           mva.ptJet2 > 30. &&
           mva.ptmiss < 40.  ) Jet2PtCuts = true;

      // #########################################
      // 2 JETS CATEGORIES
      // #########################################
      if (Jet2PtCuts) {

        bool Jet2CutsVBFPass = false;
        if (mva.mDiJet > 650. && mva.deltaEtaJets > 3.5) Jet2CutsVBFPass = true;
        
        bool Jet2CutsGFPass = false;
        if (Jet2CutsVBFPass == false &&
            mva.mDiJet > 250. && mva.ptDiMu > 50. ) Jet2CutsGFPass = true;

        if      (Jet2CutsVBFPass) {
          histsJet2CutsVBFPass  .Fill(mva,blind);
          if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
            _dimuonMassJet2CutsVBFPass          = mva.mDiMu;
            _dimuonMassJet2CutsVBFPassResSFUp   = mDiMuResSFUp;
            _dimuonMassJet2CutsVBFPassResSFDown = mDiMuResSFDown;

            if (CutMassISO) counterBDTvbfCut[40+1]++; 
            treeJet2CutsVBFPass->Fill();
          }
        }
        else if (Jet2CutsGFPass ) {
          histsJet2CutsGFPass   .Fill(mva,blind);
          if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
            _dimuonMassJet2CutsGFPass          = mva.mDiMu;
            _dimuonMassJet2CutsGFPassResSFUp   = mDiMuResSFUp;  
            _dimuonMassJet2CutsGFPassResSFDown = mDiMuResSFDown;

            if (CutMassISO) counterBDTvbfCut[40+2]++; 
            treeJet2CutsGFPass->Fill();
          }
        }

        else {
          histsJet2CutsFailVBFGF.Fill(mva,blind);
          if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
            _dimuonMassJet2CutsFailVBFGF          = mva.mDiMu;
            _dimuonMassJet2CutsFailVBFGFResSFUp   = mDiMuResSFUp;  
            _dimuonMassJet2CutsFailVBFGFResSFDown = mDiMuResSFDown;

            if (CutMassISO) counterBDTvbfCut[40+3]++; 
            treeJet2CutsFailVBFGF->Fill();
          }
        }
      } //if (Jet2PtCuts)
      

        // #########################################
        // 0+1 JET CATEGORIES
        // #########################################
      else {

        // pt(mm)> 10 GeV/c
        if (mva.ptDiMu > 10. ) {
          
          // inclusive
          histsJets01PassPtG10.Fill(mva,blind);

          if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
            _dimuonMassJets01PassPtG10          = mva.mDiMu;
            _dimuonMassJets01PassPtG10ResSFUp   = mDiMuResSFUp;  
            _dimuonMassJets01PassPtG10ResSFDown = mDiMuResSFDown;

            treeJets01PassPtG10->Fill();
          }

          // geom. categories
          if (isBB) {
            histsJets01PassPtG10BB.Fill(mva,blind);
            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01PassPtG10BB          = mva.mDiMu;
              _dimuonMassJets01PassPtG10BBResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01PassPtG10BBResSFDown = mDiMuResSFDown;

              treeJets01PassPtG10BB->Fill();
            }
          }

          if (isBO) {
            histsJets01PassPtG10BO.Fill(mva,blind);
            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01PassPtG10BO          = mva.mDiMu;
              _dimuonMassJets01PassPtG10BOResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01PassPtG10BOResSFDown = mDiMuResSFDown;

              treeJets01PassPtG10BO->Fill();
            }
          }

          if (isBE) {
            histsJets01PassPtG10BE.Fill(mva,blind);
            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01PassPtG10BE          = mva.mDiMu;
              _dimuonMassJets01PassPtG10BEResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01PassPtG10BEResSFDown = mDiMuResSFDown;

              treeJets01PassPtG10BE->Fill();
            }
          }

          if (isOO) {
            histsJets01PassPtG10OO.Fill(mva,blind);
            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01PassPtG10OO          = mva.mDiMu;
              _dimuonMassJets01PassPtG10OOResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01PassPtG10OOResSFDown = mDiMuResSFDown;

              treeJets01PassPtG10OO->Fill();
            }
          }

          if (isOE) {
            histsJets01PassPtG10OE.Fill(mva,blind);
            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01PassPtG10OE          = mva.mDiMu;
              _dimuonMassJets01PassPtG10OEResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01PassPtG10OEResSFDown = mDiMuResSFDown;

              treeJets01PassPtG10OE->Fill();
            }
          }

          if (isEE) {
            histsJets01PassPtG10EE.Fill(mva,blind);
            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01PassPtG10EE          = mva.mDiMu;
              _dimuonMassJets01PassPtG10EEResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01PassPtG10EEResSFDown = mDiMuResSFDown;

              treeJets01PassPtG10EE->Fill();
            }
          }


          if (isBE || isOO ) histsJets01PassPtG10CC.Fill(mva,blind);
          if (isOE || isEE ) histsJets01PassPtG10FF.Fill(mva,blind);
          if (CutMassISO) {
            counterDiPt10GeV[23]++;
            if (isBB) counterDiPt10GeV[23+1]++;
            if (isBO) counterDiPt10GeV[23+2]++;
            if (isBE) counterDiPt10GeV[23+3]++;
            if (isOO) counterDiPt10GeV[23+4]++;
            if (isOE) counterDiPt10GeV[23+5]++;
            if (isEE) counterDiPt10GeV[23+6]++;
            if (isBE || isOO ) counterDiPt10GeV[23+7]++;
            if (isOE || isEE ) counterDiPt10GeV[23+8]++;
          }
        }//if (mva.ptDiMu > 10. ) 

        else {
          
          // inclusive
          histsJets01FailPtG10.Fill(mva,blind);
          
          if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
            _dimuonMassJets01FailPtG10          = mva.mDiMu;
            _dimuonMassJets01FailPtG10ResSFUp   = mDiMuResSFUp;  
            _dimuonMassJets01FailPtG10ResSFDown = mDiMuResSFDown;

            treeJets01FailPtG10->Fill();
          }
          

          // geom. categories
          if (isBB) {
            histsJets01FailPtG10BB.Fill(mva,blind);

            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01FailPtG10BB          = mva.mDiMu;
              _dimuonMassJets01FailPtG10BBResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01FailPtG10BBResSFDown = mDiMuResSFDown;

              treeJets01FailPtG10BB->Fill();
            }
          }

          if (isBO) {
            histsJets01FailPtG10BO.Fill(mva,blind);

            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01FailPtG10BO          = mva.mDiMu;
              _dimuonMassJets01FailPtG10BOResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01FailPtG10BOResSFDown = mDiMuResSFDown;

              treeJets01FailPtG10BO->Fill();
            }
          }

          if (isBE) {
            histsJets01FailPtG10BE.Fill(mva,blind);

            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01FailPtG10BE          = mva.mDiMu;
              _dimuonMassJets01FailPtG10BEResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01FailPtG10BEResSFDown = mDiMuResSFDown;

              treeJets01FailPtG10BE->Fill();
            }
          }

          if (isOO) {
            histsJets01FailPtG10OO.Fill(mva,blind);

            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01FailPtG10OO          = mva.mDiMu;
              _dimuonMassJets01FailPtG10OOResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01FailPtG10OOResSFDown = mDiMuResSFDown;

              treeJets01FailPtG10OO->Fill();
            }
          }

          if (isOE) {
            histsJets01FailPtG10OE.Fill(mva,blind);

            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01FailPtG10OE          = mva.mDiMu;
              _dimuonMassJets01FailPtG10OEResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01FailPtG10OEResSFDown = mDiMuResSFDown;

              treeJets01FailPtG10OE->Fill();
            }
          }

          if (isEE) {
            histsJets01FailPtG10EE.Fill(mva,blind);

            if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {
              _dimuonMassJets01FailPtG10EE          = mva.mDiMu;
              _dimuonMassJets01FailPtG10EEResSFUp   = mDiMuResSFUp;  
              _dimuonMassJets01FailPtG10EEResSFDown = mDiMuResSFDown;

              treeJets01FailPtG10EE->Fill();
            }
          }

          if (isBE || isOO ) histsJets01FailPtG10CC.Fill(mva,blind);
          if (isOE || isEE ) histsJets01FailPtG10FF.Fill(mva,blind);

          if (CutMassISO) {
            counterDiPt10GeV[32]++;
            if (isBB) counterDiPt10GeV[32+1]++;
            if (isBO) counterDiPt10GeV[32+2]++;
            if (isBE) counterDiPt10GeV[32+3]++;
            if (isOO) counterDiPt10GeV[32+4]++;
            if (isOE) counterDiPt10GeV[32+5]++;
            if (isEE) counterDiPt10GeV[32+6]++;
            if (isBE || isOO ) counterDiPt10GeV[32+7]++;
            if (isOE || isEE ) counterDiPt10GeV[32+8]++;
          } 
        }

      }
      
    } //atLeastOneMuonFired
    
      ///////////////////////////////////////////


    if (!isHltMatched(reco1,reco2,allowedHLTPaths))
      continue;


    //4 GeV Window Plots
    if (mva.mDiMu < 127.0 && mva.mDiMu > 123.0)
      {
        if (!blind)
          {
            hists4GeVWindow.Fill(mva,blind);
          }
      }

    if (isBB)
      {
        histsBB.Fill(mva,blind);
      }

    if (isBO)
      {
        histsBO.Fill(mva,blind);
      }

    if (isBE)
      {
        histsBE.Fill(mva,blind);
      }

    if (isOO)
      {
        histsOO.Fill(mva,blind);
      }

    if (isOE)
      {
        histsOE.Fill(mva,blind);
      }

    if (isEE)
      {
        histsEE.Fill(mva,blind);
      }

    if (isNotBB)
      {
        histsNotBB.Fill(mva,blind);
      }

    //VBF Preselected Plots
    if (vbfPreselection)
      {
        histsVBFPresel.Fill(mva,blind);
      }

    if (vbfPreselection && isBB)
      {
        histsVBFPreselBB.Fill(mva,blind);
      }

    if (vbfPreselection && isNotBB)
      {
        histsVBFPreselNotBB.Fill(mva,blind);
      }

    //Inc Preselected Plots
    if (!vbfPreselection)
      {
        histsIncPresel.Fill(mva,blind);
      }

    if (!vbfPreselection && isBB)
      {
        histsIncPreselBB.Fill(mva,blind);
      }

    if (!vbfPreselection && isBO)
      {
        histsIncPreselBO.Fill(mva,blind);
      }

    if (!vbfPreselection && isBE)
      {
        histsIncPreselBE.Fill(mva,blind);
      }

    if (!vbfPreselection && isOO)
      {
        histsIncPreselOO.Fill(mva,blind);
      }

    if (!vbfPreselection && isOE)
      {
        histsIncPreselOE.Fill(mva,blind);
      }

    if (!vbfPreselection && isEE)
      {
        histsIncPreselEE.Fill(mva,blind);
      }

    if (!vbfPreselection && isNotBB)
      {
        histsIncPreselNotBB.Fill(mva,blind);
      }




    //VBF BDT Cut Plots
    if (vbfPreselection && passVBFBDTCut)
      {
        histsVBFBDTCut.Fill(mva,blind);
      }

    if (vbfPreselection && isBB && passVBFBDTCut)
      {
        histsVBFBDTCutBB.Fill(mva,blind);
      }

    if (vbfPreselection && isNotBB && passVBFBDTCut)
      {
        histsVBFBDTCutNotBB.Fill(mva,blind);
      }

    //Inc BDT Cut Plots
    if (!vbfPreselection && passIncBDTCut)
      {
        histsIncBDTCut.Fill(mva,blind);
      }

    if (!vbfPreselection && isBB && passIncBDTCut)
      {
        histsIncBDTCutBB.Fill(mva,blind);
      }

    if (!vbfPreselection && isBO && passIncBDTCut)
      {
        histsIncBDTCutBO.Fill(mva,blind);
      }

    if (!vbfPreselection && isBE && passIncBDTCut)
      {
        histsIncBDTCutBE.Fill(mva,blind);
      }

    if (!vbfPreselection && isOO && passIncBDTCut)
      {
        histsIncBDTCutOO.Fill(mva,blind);
      }

    if (!vbfPreselection && isOE && passIncBDTCut)
      {
        histsIncBDTCutOE.Fill(mva,blind);
      }

    if (!vbfPreselection && isEE && passIncBDTCut)
      {
        histsIncBDTCutEE.Fill(mva,blind);
      }

    if (!vbfPreselection && isNotBB && passIncBDTCut)
      {
        histsIncBDTCutNotBB.Fill(mva,blind);
      }


    if (!vbfPreselection && mva.ptDiMu < 20.0)
      {
        histsIncPreselDiMuPtL20.Fill(mva,blind);
      }

    if (vbfPreselection && mva.ptDiMu < 20.0)
      {
        histsVBFPreselDiMuPtL20.Fill(mva,blind);
      }

    if (vbfPreselection && mva.ptmiss < 50.0)
      {
        histsVBFPreselPtMiss50Veto.Fill(mva,blind);
      }

    if (!vbfPreselection && mva.ptDiMu > 10.0)
      {
        histsIncPreselPtG10.Fill(mva,blind);
      }

    if (!vbfPreselection && isBB && mva.ptDiMu > 10.0)
      {
        histsIncPreselPtG10BB.Fill(mva,blind);
        if(mva.RelMassRes < 0.009) {histsIncPreselPtG10BBres1.Fill(mva,blind);}
        else{histsIncPreselPtG10BBres2.Fill(mva,blind);}
        if(mva.RelMassResCov < 0.010) {histsIncPreselPtG10BBres1Cov.Fill(mva,blind);}
        else{histsIncPreselPtG10BBres2Cov.Fill(mva,blind);}
      }

    if (!vbfPreselection && isBO && mva.ptDiMu > 10.0)
      {
        histsIncPreselPtG10BO.Fill(mva,blind);
        if(mva.RelMassRes < 0.011) {histsIncPreselPtG10BOres1.Fill(mva,blind);}
        else{histsIncPreselPtG10BOres2.Fill(mva,blind);}
      }

    if (!vbfPreselection && isBE && mva.ptDiMu > 10.0)
      {
        histsIncPreselPtG10BE.Fill(mva,blind);
      }

    if (!vbfPreselection && isOO && mva.ptDiMu > 10.0)
      {
        histsIncPreselPtG10OO.Fill(mva,blind);
      }

    if (!vbfPreselection && isOE && mva.ptDiMu > 10.0)
      {
        histsIncPreselPtG10OE.Fill(mva,blind);
      }

    if (!vbfPreselection && isEE && mva.ptDiMu > 10.0)
      {
        histsIncPreselPtG10EE.Fill(mva,blind);
      }

    if (!vbfPreselection && isNotBB && mva.ptDiMu > 10.0)
      {
        histsIncPreselPtG10NotBB.Fill(mva,blind);
      }

    if (vbfPreselection && mva.mDiJet > 550.)
      {
        histsVBFMJJG550.Fill(mva,blind);
      }

    if (vbfPreselection        && 
        mva.mDiJet > 550.      && 
        mva.deltaEtaJets > 3.5 && 
        mva.ptmiss<100. )
      {
        histsVBFDeJJG3p5MJJG550pTmissL100.Fill(mva,blind);
      }

    if (vbfPreselection        && 
        mva.mDiJet > 500.      && 
        mva.deltaEtaJets > 3.4 && 
        mva.ptmiss<25. )
      {
        histsVBFDeJJG3p4MJJG500pTmissL25.Fill(mva,blind);
      }

    ///////////////////////////////////////////

    timeReading += difftime(timeStopReading,timeStartReading);
    timeProcessing += difftime(timeStartFilling,timeStopReading);
    timeFilling += difftime(time(NULL),timeStartFilling);
  }// end event loop
  time_t timeEndEventLoop = time(NULL);

  // for efficiency:
  //////////////////////////////////////////////////////////////////////
  for(unsigned iS = 0; iS<Nbin ;iS++){
    if(Flag_subCat == 1) counterGenBoson[iS] = counterGenBoson[0]; // the same dominator for sub catergories
    std::cout << " ########################################## \n";
    std::cout << " ########################################## \n";
    std::cout << " For subcategory = " << iS << " Events after: \n\n";
    //counterGenBoson = counterGenBoson - 40000;
    std::cout << " Gen Mass   Cuts = " << counterGenBoson[iS] << std::endl;
    std::cout << " Minimal    Cuts = " << counterMinimCuts[iS] << std::endl;
    std::cout << " pT/eta     Cuts = " << counterAccCuts[iS]   << std::endl;
    std::cout << " Muon ID    Cuts = " << counterIDCuts[iS]    << std::endl;
    std::cout << " Muon Iso   Cuts = " << counterIsoCuts[iS]   << std::endl;
    std::cout << " Trigger    Cuts = " << counterTrigCuts[iS]  << std::endl;
    std::cout << " VBF pres.  Cuts = " << counterJetSel[iS]  << std::endl;
    std::cout << " VBF BDT    Cuts = " << counterBDTvbfCut[iS]  << std::endl;
    std::cout << " non VBF pres.  Cuts = " << counterNonJetSel[iS]  << std::endl;
    std::cout << " pt(mumu)>10 GeV = " << counterDiPt10GeV[iS]  << std::endl;
    // calculate efficiency and binom. error:
    // calculate efficiency and binom. error:
    EffMinimCuts[iS] = float(counterMinimCuts[iS])/float(counterGenBoson[iS]);
    EffAccCuts[iS] = float(counterAccCuts[iS])/float(counterGenBoson[iS]);
    EffIDCuts[iS] = float(counterIDCuts[iS])/float(counterGenBoson[iS]);
    EffIsoCuts[iS] = float(counterIsoCuts[iS])/float(counterGenBoson[iS]);
    EffTrigCuts[iS] = float(counterTrigCuts[iS])/float(counterGenBoson[iS]);
    EffJetSel[iS] = float(counterJetSel[iS])/float(counterGenBoson[iS]);
    EffBDTvbfCut[iS] = float(counterBDTvbfCut[iS])/float(counterGenBoson[iS]);
    EffNonJetSel[iS] = float(counterNonJetSel[iS])/float(counterGenBoson[iS]);
    EffDiPt10GeV[iS] = float(counterDiPt10GeV[iS])/float(counterGenBoson[iS]);

    int iTypeDEff = 1;// 1 - TEfficiency, 2 - Benomial  
    if(iTypeDEff == 2){
      dEffMinimCuts[iS] = sqrt( EffMinimCuts[iS]*(1-EffMinimCuts[iS])/float(counterGenBoson[iS]) );
      dEffAccCuts[iS] = sqrt( EffAccCuts[iS]*(1-EffAccCuts[iS])/float(counterGenBoson[iS]) );
      dEffIDCuts[iS] = sqrt( EffIDCuts[iS]*(1-EffIDCuts[iS])/float(counterGenBoson[iS]) );
      dEffIsoCuts[iS] = sqrt( EffIsoCuts[iS]*(1-EffIsoCuts[iS])/float(counterGenBoson[iS]) );
      dEffTrigCuts[iS] = sqrt( EffTrigCuts[iS]*(1-EffTrigCuts[iS])/float(counterGenBoson[iS]) );
      dEffJetSel[iS] = sqrt( EffJetSel[iS]*(1-EffJetSel[iS])/float(counterGenBoson[iS]) );
      dEffBDTvbfCut[iS] = sqrt( EffBDTvbfCut[iS]*(1-EffBDTvbfCut[iS])/float(counterGenBoson[iS]) );
      dEffNonJetSel[iS] = sqrt( EffNonJetSel[iS]*(1-EffNonJetSel[iS])/float(counterGenBoson[iS]) );
      dEffDiPt10GeV[iS] = sqrt( EffDiPt10GeV[iS]*(1-EffDiPt10GeV[iS])/float(counterGenBoson[iS]) );
    }
    if(iTypeDEff == 1){
      //0.683 - 1 sigma, true - upper, false - lower boundary
      bool upper = true;
      dEffMinimCuts[iS] = TEfficiency::ClopperPearson( counterGenBoson[iS], counterMinimCuts[iS], 0.683, upper ) - EffMinimCuts[iS];
      dEffAccCuts[iS]   = TEfficiency::ClopperPearson( counterGenBoson[iS], counterAccCuts[iS],   0.683, upper ) - EffAccCuts[iS];
      dEffIDCuts[iS]    = TEfficiency::ClopperPearson( counterGenBoson[iS], counterIDCuts[iS],    0.683, upper ) - EffIDCuts[iS];
      dEffIsoCuts[iS]   = TEfficiency::ClopperPearson( counterGenBoson[iS], counterIsoCuts[iS],   0.683, upper ) - EffIsoCuts[iS];
      dEffTrigCuts[iS]  = TEfficiency::ClopperPearson( counterGenBoson[iS], counterTrigCuts[iS],  0.683, upper ) - EffTrigCuts[iS];
      dEffJetSel[iS]    = TEfficiency::ClopperPearson( counterGenBoson[iS], counterJetSel[iS],    0.683, upper ) - EffJetSel[iS];
      dEffBDTvbfCut[iS] = TEfficiency::ClopperPearson( counterGenBoson[iS], counterBDTvbfCut[iS], 0.683, upper ) - EffBDTvbfCut[iS];
      dEffNonJetSel[iS] = TEfficiency::ClopperPearson( counterGenBoson[iS], counterNonJetSel[iS], 0.683, upper ) - EffNonJetSel[iS];
      dEffDiPt10GeV[iS] = TEfficiency::ClopperPearson( counterGenBoson[iS], counterDiPt10GeV[iS], 0.683, upper ) - EffDiPt10GeV[iS];
    }
    std::cout << " ########################################## \n";
    std::cout << " Efficiency after selection: \n\n";
    cout.unsetf(ios::floatfield);            // floatfield not set
    cout.precision(3);
    std::cout << " Opposit charge, M_RECO = 110-150 GeV  = " << EffMinimCuts[iS] <<" +/- ";
    cout.precision(1);
    std::cout << dEffMinimCuts[iS] << std::endl;
    cout.precision(3);
    std::cout << " + pT > 25 GeV,|eta| < 2.1         Cut = " << EffAccCuts[iS]   <<" +/- ";
    cout.precision(1);
    std::cout << dEffAccCuts[iS]   << std::endl;
    cout.precision(3);
    std::cout << " + Tight Muon ID                   Cut = " << EffIDCuts[iS]    <<" +/- " ;
    cout.precision(1);
    std::cout <<  dEffIDCuts[iS]    << std::endl;
    cout.precision(3);
    std::cout << " + Muon Relative PF Isolation      Cut = " << EffIsoCuts[iS]   <<" +/- " ;
    cout.precision(1);
    std::cout << dEffIsoCuts[iS]   << std::endl;
    cout.precision(3);
    std::cout << " + Trigger HLT_Mu24Iso_eta2p1      Cut = " << EffTrigCuts[iS]  <<" +/- " ;
    cout.precision(1);
    std::cout << dEffTrigCuts[iS]  << std::endl;
    cout.precision(3);
    std::cout << " + VBF Jet preselection            Cut = " << EffJetSel[iS]    <<" +/- " ;
    cout.precision(1);
    std::cout << dEffJetSel[iS] << std::endl;
    cout.precision(3);
    std::cout << " + VBF BDT                         Cut = " << EffBDTvbfCut[iS]    <<" +/- " ;
    cout.precision(1);
    std::cout << dEffBDTvbfCut[iS] << std::endl;
    cout.precision(3);
    std::cout << " + non VBF Jet preselection         Cut = " << EffNonJetSel[iS]    <<" +/- " ;
    cout.precision(1);
    std::cout << dEffNonJetSel[iS] << std::endl;
    cout.precision(3);
    std::cout << " + pt(mumu)> 10 GeV/c, no VBF pres. and no BDT  Cut = " << EffDiPt10GeV[iS]    <<" +/- " ;
    cout.precision(1);
    std::cout << dEffDiPt10GeV[iS] << std::endl;

  } // end for 


  // write BDT cut efficiency to myfile
  myfile << "EffBDT{}  = {";
  myfile.precision(3);
  for(unsigned iS = 0; iS<Nbin ;iS++){
    myfile << EffBDTvbfCut[iS] << ", ";
  }
  myfile << "};\n\n";

  myfile << "dEffBDT{} = {";
  myfile.precision(1);
  for(unsigned iS = 0; iS<Nbin ;iS++){
    myfile << dEffBDTvbfCut[iS] << ", ";
  }
  myfile << "};\n\n";


  myfile.close();

  std::cout << "write efficiency for sub categories\n" << std::fixed;
  // write efficiency for sub categories to myfileSubCat
  for(unsigned iS = 0; iS < 6 ;iS++){
    myfileSubCat.precision(8);
    myfileSubCat << subCategory[iS] << std::fixed << EffDiPt10GeV[iS];
    myfileSubCat.precision(8);
    myfileSubCat << "   " << std::fixed <<  dEffDiPt10GeV[iS] << "\n";
  }
  myfileSubCat.precision(8);
  myfileSubCat << subCategory[6] << std::fixed << EffNonJetSel[6];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffNonJetSel[6] << "\n";

  myfileSubCat.precision(8);
  myfileSubCat << subCategory[7] <<  std::fixed << EffDiPt10GeV[7];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffDiPt10GeV[7] << "\n";

  myfileSubCat.precision(8);
  myfileSubCat << subCategory[8] <<  std::fixed << EffJetSel[8];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffJetSel[8] << "\n";

  myfileSubCat.precision(8);
  myfileSubCat << subCategory[9] <<  std::fixed << EffBDTvbfCut[9];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffBDTvbfCut[9] << "\n";

  for(unsigned iS = 0; iS < 6 ;iS++){
    myfileSubCat.precision(8);
    myfileSubCat << subCategory[iS+10] <<  std::fixed << EffNonJetSel[iS];
    myfileSubCat.precision(8);
    myfileSubCat << "   " <<  std::fixed << dEffNonJetSel[iS] << "\n";
  }

  myfileSubCat.precision(8);
  myfileSubCat << subCategory[16] <<  std::fixed << EffJetSel[9];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffJetSel[9] << "\n";

  // resolution categories for BB and BO

  myfileSubCat.precision(8);
  myfileSubCat << subCategory[17] <<  std::fixed << EffDiPt10GeV[8];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffDiPt10GeV[8] << "\n";
 
  myfileSubCat.precision(8);
  myfileSubCat << subCategory[18] <<  std::fixed << EffDiPt10GeV[9];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffDiPt10GeV[9] << "\n";
 
  myfileSubCat.precision(8);
  myfileSubCat << subCategory[19] <<  std::fixed << EffDiPt10GeV[10];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffDiPt10GeV[10] << "\n";
 
  myfileSubCat.precision(8);
  myfileSubCat << subCategory[20] <<  std::fixed << EffDiPt10GeV[11];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffDiPt10GeV[11] << "\n";
 
  // VBF cutbased efficiency

  myfileSubCat.precision(8);
  myfileSubCat << subCategory[21] <<  std::fixed << EffJetSel[10];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffJetSel[10] << "\n";

  myfileSubCat.precision(8);
  myfileSubCat << subCategory[22] <<  std::fixed << EffJetSel[11];
  myfileSubCat.precision(8);
  myfileSubCat << "   " <<  std::fixed << dEffJetSel[11] << "\n";

  //Baseline++

  // write efficiency for sub categories to myfileSubCat
  for(unsigned iS = 23; iS < 41 ;iS++){
    myfileSubCat.precision(8);
    myfileSubCat << subCategory[iS] <<  std::fixed << EffDiPt10GeV[iS];
    myfileSubCat.precision(8);
    myfileSubCat << "   " <<  std::fixed << dEffDiPt10GeV[iS] << "\n";
  }
 
  for(unsigned iS = 41; iS < 44 ;iS++){
    myfileSubCat.precision(8);
    myfileSubCat << subCategory[iS] <<  std::fixed << EffBDTvbfCut[iS];
    myfileSubCat.precision(8);
    myfileSubCat << "   " <<  std::fixed << dEffBDTvbfCut[iS] << "\n";
  }
 

  myfileSubCat.close();


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  // write the ttree
  cout << "About to begin writting stuff to files..." << endl;
  outFile->cd();
  cout << "Writing tree..." << endl;
  _outTree->Write();

  cout << "Writing first hists..." << endl;
  hists.Write(outFile,"");
  /*
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

    histsIncBDTCut.Write(outFile,"IncBDTCut");
    histsIncBDTCutBB.Write(outFile,"IncBDTCutBB");
    histsIncBDTCutBO.Write(outFile,"IncBDTCutBO");
    histsIncBDTCutBE.Write(outFile,"IncBDTCutBE");
    histsIncBDTCutOO.Write(outFile,"IncBDTCutOO");
    histsIncBDTCutOE.Write(outFile,"IncBDTCutOE");
    histsIncBDTCutEE.Write(outFile,"IncBDTCutEE");
    histsIncBDTCutNotBB.Write(outFile,"IncBDTCutNotBB");

    histsVBFBDTCut.Write(outFile,"VBFBDTCut");
    histsVBFBDTCutBB.Write(outFile,"VBFBDTCutBB");
    histsVBFBDTCut.Write(outFile,"VBFBDTCutNotBB");

    histsVBFPreselDiMuPtL20.Write(outFile,"VBFPreselDiMuPtL20");
    histsIncPreselDiMuPtL20.Write(outFile,"IncPreselDiMuPtL20");

    histsIncPreselPUJETID.Write(outFile,"IncPreselPUJETID");
    histsVBFPreselPUJETID.Write(outFile,"VBFPreselPUJETID");

    histsIncPreselPUJETIDForVeto.Write(outFile,"IncPreselPUJETIDForVeto");
    histsVBFPreselPUJETIDForVeto.Write(outFile,"VBFPreselPUJETIDForVeto");

    histsVBFPreselPtMiss50Veto.Write(outFile,"VBFPreselPtMiss50Veto");

    histsIncPreselPtG10.Write(outFile,"IncPreselPtG10");
    histsIncPreselPtG10BB.Write(outFile,"IncPreselPtG10BB");
    histsIncPreselPtG10BO.Write(outFile,"IncPreselPtG10BO");
    histsIncPreselPtG10BE.Write(outFile,"IncPreselPtG10BE");
    histsIncPreselPtG10OO.Write(outFile,"IncPreselPtG10OO");
    histsIncPreselPtG10OE.Write(outFile,"IncPreselPtG10OE");
    histsIncPreselPtG10EE.Write(outFile,"IncPreselPtG10EE");
    histsIncPreselPtG10NotBB.Write(outFile,"IncPreselPtG10NotBB");
    histsIncPreselPtG10BBres1.Write(outFile,"IncPreselPtG10BBres1");
    histsIncPreselPtG10BBres2.Write(outFile,"IncPreselPtG10BBres2");
    histsIncPreselPtG10BOres1.Write(outFile,"IncPreselPtG10BOres1");
    histsIncPreselPtG10BOres2.Write(outFile,"IncPreselPtG10BOres2");
    histsIncPreselPtG10BBres1Cov.Write(outFile,"IncPreselPtG10BBres1Cov");
    histsIncPreselPtG10BBres2Cov.Write(outFile,"IncPreselPtG10BBres2Cov");

    histsVBFMJJG550.Write(outFile,"VBFMJJG550");
    histsVBFDeJJG3p5MJJG550pTmissL100.Write(outFile,"histsVBFDeJJG3p5MJJG550pTmissL100");
    histsVBFDeJJG3p4MJJG500pTmissL25 .Write(outFile,"histsVBFDeJJG3p4MJJG500pTmissL25");
  */

  cout << "Writing baseline++..." << endl;
  // baseline++
  histsJets01PassPtG10  .Write(outFile,"Jets01PassPtG10");
  histsJets01PassPtG10BB.Write(outFile,"Jets01PassPtG10BB");
  histsJets01PassPtG10BO.Write(outFile,"Jets01PassPtG10BO");
  histsJets01PassPtG10BE.Write(outFile,"Jets01PassPtG10BE");
  histsJets01PassPtG10OO.Write(outFile,"Jets01PassPtG10OO");
  histsJets01PassPtG10OE.Write(outFile,"Jets01PassPtG10OE");
  histsJets01PassPtG10EE.Write(outFile,"Jets01PassPtG10EE");
  histsJets01PassPtG10CC.Write(outFile,"Jets01PassPtG10CC");
  histsJets01PassPtG10FF.Write(outFile,"Jets01PassPtG10FF");

  histsJets01FailPtG10  .Write(outFile,"Jets01FailPtG10");
  histsJets01FailPtG10BB.Write(outFile,"Jets01FailPtG10BB");
  histsJets01FailPtG10BO.Write(outFile,"Jets01FailPtG10BO");
  histsJets01FailPtG10BE.Write(outFile,"Jets01FailPtG10BE");
  histsJets01FailPtG10OO.Write(outFile,"Jets01FailPtG10OO");
  histsJets01FailPtG10OE.Write(outFile,"Jets01FailPtG10OE");
  histsJets01FailPtG10EE.Write(outFile,"Jets01FailPtG10EE");
  histsJets01FailPtG10CC.Write(outFile,"Jets01FailPtG10CC");
  histsJets01FailPtG10FF.Write(outFile,"Jets01FailPtG10FF");

  histsJet2CutsVBFPass.Write(outFile,"Jet2CutsVBFPass");
  histsJet2CutsGFPass.Write(outFile,"Jet2CutsGFPass");
  histsJet2CutsFailVBFGF.Write(outFile,"Jet2CutsFailVBFGF");

  treeJets01PassPtG10   -> Write();
  treeJets01PassPtG10BB -> Write();
  treeJets01PassPtG10BO -> Write();
  treeJets01PassPtG10BE -> Write();
  treeJets01PassPtG10OO -> Write();
  treeJets01PassPtG10OE -> Write();
  treeJets01PassPtG10EE -> Write();

  treeJets01FailPtG10   -> Write();
  treeJets01FailPtG10BB -> Write();
  treeJets01FailPtG10BO -> Write();
  treeJets01FailPtG10BE -> Write();
  treeJets01FailPtG10OO -> Write();
  treeJets01FailPtG10OE -> Write();
  treeJets01FailPtG10EE -> Write();

  treeJet2CutsVBFPass   -> Write();
  treeJet2CutsGFPass    -> Write();
  treeJet2CutsFailVBFGF -> Write();

  


  //
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

  if(getPFRelIso(mu1) > 0.12 || getPFRelIso(mu2) > 0.12)
    {
      //cout << "Iso 1: "<< getPFRelIso(mu1) << "    Iso 2: " << getPFRelIso(mu2) << endl;
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

  if(getPFRelIso(mu1) > 0.12 || getPFRelIso(mu2) > 0.12)
    {
      //cout << "Iso 1: "<< getPFRelIso(mu1) << "    Iso 2: " << getPFRelIso(mu2) << endl;
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

  unsigned nMassBins = 800;
  float minMass = 0.;
  float maxMass = 400.;
  unsigned nMVABins = 200;

  mDiMu = new TH1F("mDiMu","DiMuon Mass",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMu);

  mDiMu110to160 = new TH1F("mDiMu110to160","DiMuon Mass in [110,160]",100,110.,160.);
  histVec.push_back(mDiMu110to160);

  mDiMuTrkRelIso = new TH1F("mDiMuTrkRelIso","DiMuon Mass w/ TrkRelIso",nMassBins,minMass,maxMass);
  histVec.push_back(mDiMuTrkRelIso);

  RelMassRes = new TH1F("RelMassRes","Mass Resoluton",100,0.,0.05);
  histVec.push_back(RelMassRes);
  RelMassResCov = new TH1F("RelMassResCov","Mass Resoluton",100,0.,0.05);
  histVec.push_back(RelMassResCov);

  mDiJet = new TH1F("mDiJet","DiJet Mass",500,0,2000);
  histVec.push_back(mDiJet);

  ptDiMu = new TH1F("ptDiMu","DiMuon Pt",250,0,500);
  histVec.push_back(ptDiMu);

  ptDiJet = new TH1F("ptDiJet","DiJet Pt",250,0,1000);
  histVec.push_back(ptDiJet);

  yDiMu = new TH1F("yDiMu","DiMuon Rapidity",100,-4,4);
  histVec.push_back(yDiMu);
  yDiJet = new TH1F("yDiJet","DiJet Rapidity",100,-5,5);
  histVec.push_back(yDiJet);

  yVptDiMu = new TH2F("yVptDiMu","DiMuon Rapidity v. p_{T}",250,0,500,100,0,4);
  histVec2D.push_back(yVptDiMu);
  ptVmDiMu = new TH2F("ptVmDiMu","DiMuon p_{T} v. Mass",nMassBins,minMass,maxMass,250,0,250);
  histVec2D.push_back(ptVmDiMu);
  yVmDiMu = new TH2F("yVmDiMu","DiMuon |y| v. Mass",nMassBins,minMass,maxMass,100,0,4);
  histVec2D.push_back(yVmDiMu);
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
  deltaPhiHJ1 = new TH1F("deltaPhiHJ1","#Delta #phi Leading Jet Dimuon",50,0.0,3.2);
  histVec.push_back(deltaPhiHJ1);

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

  BDTHistVBF = new TH1F("BDTHistVBF","BDT Discriminator",nMVABins,-1,1);
  histVec.push_back(BDTHistVBF);

  BDTHistMuonOnlyVMass = new TH2F("BDTHistMuonOnlyVMass","BDT Discriminator",nMassBins,minMass,maxMass,nMVABins,-1,1);
  histVec2D.push_back(BDTHistMuonOnlyVMass);
  BDTHistVBFVMass = new TH2F("BDTHistVBFVMass","BDT Discriminator",nMassBins,minMass,maxMass,nMVABins,-1,1);
  histVec2D.push_back(BDTHistVBFVMass);

  relIsoMu1 = new TH1F("relIsoMu1","",1000,0,10.0);
  histVec.push_back(relIsoMu1);
  relIsoMu2 = new TH1F("relIsoMu2","",1000,0,10.0);
  histVec.push_back(relIsoMu2);

  trkRelIsoMu1 = new TH1F("trkRelIsoMu1","",1000,0,10.0);
  histVec.push_back(trkRelIsoMu1);
  trkRelIsoMu2 = new TH1F("trkRelIsoMu2","",1000,0,10.0);
  histVec.push_back(trkRelIsoMu2);

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
  ptmiss = new TH1F("ptmiss","",160,0,800);
  histVec.push_back(ptmiss);
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

void 
HistStruct::Fill(const MVA& mva, bool blind)
{

  if (mva.trkRelIsoMu1 < 0.1 && mva.trkRelIsoMu2 < 0.1) {
    trkRelIsoMu1->Fill(mva.trkRelIsoMu1, mva.weight);
    trkRelIsoMu2->Fill(mva.trkRelIsoMu2, mva.weight);
 
    if (!blind) mDiMuTrkRelIso->Fill(mva.mDiMu, mva.weight);
  } // is Trk Rel Iso passed on both muons

    
  if (mva.relIsoMu1 < 0.12 && mva.relIsoMu2 < 0.12) {

    yDiMu->Fill(mva.yDiMu, mva.weight);
    ptDiMu->Fill(mva.ptDiMu, mva.weight);
    ptMu1->Fill(mva.ptMu1, mva.weight);
    ptMu2->Fill(mva.ptMu2, mva.weight);
    etaMu1->Fill(mva.etaMu1, mva.weight);
    etaMu2->Fill(mva.etaMu2, mva.weight);
    cosThetaStar->Fill(mva.cosThetaStar, mva.weight);
    cosThetaStarCS->Fill(mva.cosThetaStarCS, mva.weight);
    deltaPhiMuons->Fill(mva.deltaPhiMuons, mva.weight);
    deltaEtaMuons->Fill(mva.deltaEtaMuons, mva.weight);
    deltaRMuons->Fill(mva.deltaRMuons, mva.weight);
    relIsoMu1->Fill(mva.relIsoMu1, mva.weight);
    relIsoMu2->Fill(mva.relIsoMu2, mva.weight);
    nPU->Fill(mva.nPU, mva.weight);
    nVtx->Fill(mva.nVtx, mva.weight);
    met->Fill(mva.met, mva.weight);
    ptmiss->Fill(mva.ptmiss, mva.weight);
    
    mDiJet->Fill(mva.mDiJet, mva.weight);
    ptDiJet->Fill(mva.ptDiJet, mva.weight);
    yDiJet->Fill(mva.yDiJet, mva.weight);
    ptJet1->Fill(mva.ptJet1, mva.weight);
    ptJet2->Fill(mva.ptJet2, mva.weight);
    etaJet1->Fill(mva.etaJet1, mva.weight);
    etaJet2->Fill(mva.etaJet2, mva.weight);
    deltaEtaJets->Fill(mva.deltaEtaJets, mva.weight);
    deltaPhiJets->Fill(mva.deltaPhiJets, mva.weight);
    deltaRJets->Fill(mva.deltaRJets, mva.weight);
    deltaPhiHJ1->Fill(mva.deltaPhiHJ1, mva.weight);
    nJetsInRapidityGap->Fill(mva.nJetsInRapidityGap, mva.weight);
    htInRapidityGap->Fill(mva.htInRapidityGap, mva.weight);
    nJets->Fill(mva.nJets, mva.weight);
    ht->Fill(mva.ht, mva.weight);
    
    puJetIDSimpleDiscJet1->Fill(mva.puJetIDSimpleDiscJet1,mva.weight);
    puJetIDSimpleDiscJet2->Fill(mva.puJetIDSimpleDiscJet2,mva.weight);
    puJetIDSimpleDiscJet3->Fill(mva.puJetIDSimpleDiscJet3,mva.weight);
    
    puJetIDSimpleJet1->Fill(mva.puJetIDSimpleJet1,mva.weight);
    puJetIDSimpleJet2->Fill(mva.puJetIDSimpleJet2,mva.weight);
    puJetIDSimpleJet3->Fill(mva.puJetIDSimpleJet3,mva.weight);
    
    yVptDiMu->Fill(mva.ptDiMu,fabs(mva.yDiMu),mva.weight);
    
    if (!blind)
      {

        mDiMu->Fill(mva.mDiMu, mva.weight);

        if (mva.mDiMu >= 110. && mva.mDiMu <= 160.)
          mDiMu110to160 ->Fill(mva.mDiMu, mva.weight);

        RelMassRes->Fill(mva.RelMassRes, mva.weight);
        RelMassResCov->Fill(mva.RelMassResCov, mva.weight);

        yVmDiMu->Fill(mva.mDiMu,fabs(mva.yDiMu),mva.weight);
        ptVmDiMu->Fill(mva.mDiMu,mva.ptDiMu,mva.weight);
        
        if(mva.nJets < 2)
          {
            BDTHistMuonOnly->Fill(mva.bdtValInc, mva.weight);
            BDTHistMuonOnlyVMass->Fill(mva.mDiMu, mva.bdtValInc, mva.weight);
          }
        else
          {
            BDTHistVBF->Fill(mva.bdtValVBF, mva.weight);
            BDTHistVBFVMass->Fill(mva.mDiMu, mva.bdtValVBF, mva.weight);
          }
      }

  } // is PF Rel Iso passed on both muons
  
}

double GetMassRes(_MuonInfo& mu1, _MuonInfo& mu2) {

  //double const MASS_MUON = 0.105658367;    //GeV/c2

  // get the dimuon candidate
  //TLorentzVector dimuon = GetLorentzVector(pair);
  TLorentzVector mu1Vec;
  TLorentzVector mu2Vec;
  mu1Vec.SetPtEtaPhiM(mu1.pt,mu1.eta,mu1.phi,MASS_MUON);
  mu2Vec.SetPtEtaPhiM(mu2.pt,mu2.eta,mu2.phi,MASS_MUON);
  TLorentzVector dimuon = mu1Vec + mu2Vec;

  // get the dimuon mass
  
  double dimuonMass = dimuon.M();


  TLorentzVector muon1withError, muon2, dimuon_var1; 
  muon1withError.SetPtEtaPhiM(mu1.pt+mu1.ptErr, 
                              mu1.eta, mu1.phi, MASS_MUON);
  muon2.SetPtEtaPhiM(mu2.pt, mu2.eta, mu2.phi, MASS_MUON);
  dimuon_var1 = muon1withError+muon2;
  double deltaM1 = fabs(dimuon_var1.M() - dimuonMass );

  TLorentzVector muon1, muon2withError, dimuon_var2; 
  muon1.SetPtEtaPhiM(mu1.pt, mu1.eta, mu1.phi, MASS_MUON);
  muon2withError.SetPtEtaPhiM(mu2.pt+mu2.ptErr, 
                              mu2.eta, mu2.phi, MASS_MUON);
  dimuon_var2 = muon1+muon2withError;
  double deltaM2 = fabs(dimuon_var2.M() - dimuonMass );

  return sqrt( (deltaM1*deltaM1) + (deltaM2*deltaM2) );

}
