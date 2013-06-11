#include <iomanip>
#include <map>
#include <set>
#include <vector>
#include <string>
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

#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include "boost/algorithm/string.hpp"
#include <boost/lexical_cast.hpp>

#include <limits.h>

#include <boost/lexical_cast.hpp>
//Defines method of std::string that appends any type :-)
#define appendAny(a) append(boost::lexical_cast<std::string>(a))

using namespace std;
using namespace boost;

int main(int argc, char *argv[])
{
  cout.precision(5);
  /////////////////////////////////////////////
  //////////// Configuration Options //////////
  /////////////////////////////////////////////

  const char* optionIntro = "H->MuMu Analyzer\n\nUsage: ./analyzer [--help] [--train] [--maxEvents N] <outputFileName.root> <inputFileName.root> [<inputFileName2.root>...]\n\nAllowed Options";
  program_options::options_description optionDesc(optionIntro);
  optionDesc.add_options()
      ("help,h", "Produce Help Message")
      ("outfilename",program_options::value<std::vector<string> >(), "Output File Name")
      ("firstData,f",program_options::value<vector<string> >(), "Filenames for First Dataset")
      ("secondData,s",program_options::value<vector<string> >(), "Filenames for Second Dataset")
      ("useMuscleFitSecond", "Uses Musclefit for cuts on second dataset.  Normally, use musclefit for cuts on first dataset, but not for second.")
      ("trkIso,t","Select with trk iso instead of PF iso")
  ;
  
  program_options::positional_options_description optionPos;
  optionPos.add("outfilename",-1);
  
  program_options::variables_map optionMap;
  program_options::store(program_options::command_line_parser(argc, argv).options(optionDesc).positional(optionPos).run(), optionMap);
  program_options::notify(optionMap);    
  
  if (optionMap.count("help")) 
  {
      cout << optionDesc << "\n";
      return 1;
  }

  vector<string> firstDataFileNames;
  vector<string>::const_iterator filenameItr;
  if (optionMap.count("firstData")>0)
  {
     firstDataFileNames = optionMap["firstData"].as<vector<string> >();
     if(firstDataFileNames.size()<1)
     {
       cout << "Error: Need first dataset input file names, exiting." << endl;
       return 1;
     }
  }
  else
  {
     cout << "Error: First dataset input file name arguments required, exiting." << endl;
     return 1;
  }

  vector<string> secondDataFileNames;
  if (optionMap.count("secondData")>0)
  {
     secondDataFileNames = optionMap["secondData"].as<vector<string> >();
     if(secondDataFileNames.size()<1)
     {
       cout << "Error: Need second dataset input file names, exiting." << endl;
       return 1;
     }
  }
  else
  {
     cout << "Error: Second dataset input file name arguments required, exiting." << endl;
     return 1;
  }

  string outFileName;
  if (optionMap.count("outfilename")>0)
  {
     vector<string> outFileNames = optionMap["outfilename"].as<vector<string> >();

     if (outFileNames.size() != 1)
     {
       cout << "Error: one outfilename  required, exiting." << endl;
       return 1;
     }
     outFileName = outFileNames[0];
  }
  else
  {
     cout << "Error: one outfilename  required, exiting." << endl;
     return 1;
  }

  bool useMuscleFitSecond=false;
  if (optionMap.count("useMuscleFitSecond")) 
  {
    useMuscleFitSecond=true;
    cout << "Using Musclefit for cuts on second dataset";
  }

  bool useTrkIso=false;
  if (optionMap.count("trkIso")) 
  {
    useTrkIso=true;
  }

  TChain* tree1 = new TChain("outtree");
  TChain* tree2 = new TChain("outtree");

  for(filenameItr = firstDataFileNames.begin();filenameItr != firstDataFileNames.end();filenameItr++)
  {
    tree1->AddFile(filenameItr->c_str());
  }

  for(filenameItr = secondDataFileNames.begin();filenameItr != secondDataFileNames.end();filenameItr++)
  {
    tree2->AddFile(filenameItr->c_str());
  }

  TFile* outfile = new TFile(outFileName.c_str(),"RECREATE");
  TTree* matchTree = new TTree("match","");
  TTree* only1Tree = new TTree("only1","");
  TTree* only2Tree = new TTree("only2","");
  TTree* all1Tree = new TTree("all1","");
  TTree* all2Tree = new TTree("all2","");

  map<string,int> events1;
  map<string,int> events2;
  map<string,int>::const_iterator events1Itr;
  map<string,int>::const_iterator events2Itr;

  set<int> runs1;
  set<int> runs2;
  set<int>::const_iterator runs1Itr;
  set<int>::const_iterator runs2Itr;

  int eventInfo_run1;
  int eventInfo_event1;
  TBranch* eventInfo_run1Branch = tree1->GetBranch("eventInfo_run");
  TBranch* eventInfo_event1Branch = tree1->GetBranch("eventInfo_event");
  tree1->SetBranchAddress("eventInfo_run", &eventInfo_run1);
  tree1->SetBranchAddress("eventInfo_event", &eventInfo_event1);

  float dimuonMass1;
  float dimuonPt1;
  float dimuonMass_noMuscle1;
  float dimuonPt_noMuscle1;
  tree1->SetBranchAddress("dimuonMass", &dimuonMass1);
  tree1->SetBranchAddress("dimuonPt", &dimuonPt1);
  tree1->SetBranchAddress("dimuonMass_noMuscle", &dimuonMass_noMuscle1);
  tree1->SetBranchAddress("dimuonPt_noMuscle", &dimuonPt_noMuscle1);

  int eventInfo_run2;
  int eventInfo_event2;
  TBranch* eventInfo_run2Branch = tree2->GetBranch("eventInfo_run");
  TBranch* eventInfo_event2Branch = tree2->GetBranch("eventInfo_event");
  tree2->SetBranchAddress("eventInfo_run", &eventInfo_run2);
  tree2->SetBranchAddress("eventInfo_event", &eventInfo_event2);

  float dimuonMass2;
  float dimuonPt2;
  float dimuonMass_noMuscle2;
  float dimuonPt_noMuscle2;
  tree2->SetBranchAddress("dimuonMass", &dimuonMass2);
  tree2->SetBranchAddress("dimuonPt", &dimuonPt2);
  tree2->SetBranchAddress("dimuonMass_noMuscle", &dimuonMass_noMuscle2);
  tree2->SetBranchAddress("dimuonPt_noMuscle", &dimuonPt_noMuscle2);

  float muonLead_pt1;
  float muonLead_pt_noMuscle1;
  float muonLead_eta1;
  float muonLead_phi1;
  int muonLead_passTrkRelIso1;
  int muonLead_passPFRelIso1;
  int muonLead_passTrkRelIso_noMuscle1;
  int muonLead_passPFRelIso_noMuscle1;
  int muonLead_isHltMatched1;
  tree1->SetBranchAddress("muonLead_pt", &muonLead_pt1);
  tree1->SetBranchAddress("muonLead_pt_noMuscle", &muonLead_pt_noMuscle1);
  tree1->SetBranchAddress("muonLead_eta", &muonLead_eta1);
  tree1->SetBranchAddress("muonLead_phi", &muonLead_phi1);
  tree1->SetBranchAddress("muonLead_passTrkRelIso", &muonLead_passTrkRelIso1);
  tree1->SetBranchAddress("muonLead_passPFRelIso", &muonLead_passPFRelIso1);
  tree1->SetBranchAddress("muonLead_passTrkRelIso_noMuscle", &muonLead_passTrkRelIso_noMuscle1);
  tree1->SetBranchAddress("muonLead_passPFRelIso_noMuscle", &muonLead_passPFRelIso_noMuscle1);
  tree1->SetBranchAddress("muonLead_isHltMatched", &muonLead_isHltMatched1);
  float muonSub_pt1;
  float muonSub_pt_noMuscle1;
  float muonSub_eta1;
  float muonSub_phi1;
  int muonSub_passTrkRelIso1;
  int muonSub_passPFRelIso1;
  int muonSub_passTrkRelIso_noMuscle1;
  int muonSub_passPFRelIso_noMuscle1;
  int muonSub_isHltMatched1;
  tree1->SetBranchAddress("muonSub_pt", &muonSub_pt1);
  tree1->SetBranchAddress("muonSub_pt_noMuscle", &muonSub_pt_noMuscle1);
  tree1->SetBranchAddress("muonSub_eta", &muonSub_eta1);
  tree1->SetBranchAddress("muonSub_phi", &muonSub_phi1);
  tree1->SetBranchAddress("muonSub_passTrkRelIso", &muonSub_passTrkRelIso1);
  tree1->SetBranchAddress("muonSub_passPFRelIso", &muonSub_passPFRelIso1);
  tree1->SetBranchAddress("muonSub_passTrkRelIso_noMuscle", &muonSub_passTrkRelIso_noMuscle1);
  tree1->SetBranchAddress("muonSub_passPFRelIso_noMuscle", &muonSub_passPFRelIso_noMuscle1);
  tree1->SetBranchAddress("muonSub_isHltMatched", &muonSub_isHltMatched1);

  float muonLead_pt2;
  float muonLead_pt_noMuscle2;
  float muonLead_eta2;
  float muonLead_phi2;
  int muonLead_passTrkRelIso2;
  int muonLead_passPFRelIso2;
  int muonLead_passTrkRelIso_noMuscle2;
  int muonLead_passPFRelIso_noMuscle2;
  int muonLead_isHltMatched2;
  tree2->SetBranchAddress("muonLead_pt", &muonLead_pt2);
  tree2->SetBranchAddress("muonLead_pt_noMuscle", &muonLead_pt_noMuscle2);
  tree2->SetBranchAddress("muonLead_eta", &muonLead_eta2);
  tree2->SetBranchAddress("muonLead_phi", &muonLead_phi2);
  tree2->SetBranchAddress("muonLead_passTrkRelIso", &muonLead_passTrkRelIso2);
  tree2->SetBranchAddress("muonLead_passPFRelIso", &muonLead_passPFRelIso2);
  tree2->SetBranchAddress("muonLead_passTrkRelIso_noMuscle", &muonLead_passTrkRelIso_noMuscle2);
  tree2->SetBranchAddress("muonLead_passPFRelIso_noMuscle", &muonLead_passPFRelIso_noMuscle2);
  tree2->SetBranchAddress("muonLead_isHltMatched", &muonLead_isHltMatched2);
  float muonSub_pt2;
  float muonSub_pt_noMuscle2;
  float muonSub_eta2;
  float muonSub_phi2;
  int muonSub_passTrkRelIso2;
  int muonSub_passPFRelIso2;
  int muonSub_passTrkRelIso_noMuscle2;
  int muonSub_passPFRelIso_noMuscle2;
  int muonSub_isHltMatched2;
  tree2->SetBranchAddress("muonSub_pt", &muonSub_pt2);
  tree2->SetBranchAddress("muonSub_pt_noMuscle", &muonSub_pt_noMuscle2);
  tree2->SetBranchAddress("muonSub_eta", &muonSub_eta2);
  tree2->SetBranchAddress("muonSub_phi", &muonSub_phi2);
  tree2->SetBranchAddress("muonSub_passTrkRelIso", &muonSub_passTrkRelIso2);
  tree2->SetBranchAddress("muonSub_passPFRelIso", &muonSub_passPFRelIso2);
  tree2->SetBranchAddress("muonSub_passTrkRelIso_noMuscle", &muonSub_passTrkRelIso_noMuscle2);
  tree2->SetBranchAddress("muonSub_passPFRelIso_noMuscle", &muonSub_passPFRelIso_noMuscle2);
  tree2->SetBranchAddress("muonSub_isHltMatched", &muonSub_isHltMatched2);

  int eventType1;
  tree1->SetBranchAddress("eventType", &eventType1);

  int eventType2;
  tree2->SetBranchAddress("eventType", &eventType2);

  // Now Out Trees;

  matchTree->Branch("dimuonMass1",&dimuonMass1,"dimuonMass1/F");
  only1Tree->Branch("dimuonMass1",&dimuonMass1,"dimuonMass1/F");
  all1Tree->Branch("dimuonMass1",&dimuonMass1,"dimuonMass1/F");
  matchTree->Branch("dimuonMass2",&dimuonMass2,"dimuonMass2/F");
  only2Tree->Branch("dimuonMass2",&dimuonMass2,"dimuonMass2/F");
  all2Tree->Branch("dimuonMass2",&dimuonMass2,"dimuonMass2/F");

  matchTree->Branch("dimuonPt1",&dimuonPt1,"dimuonPt1/F");
  only1Tree->Branch("dimuonPt1",&dimuonPt1,"dimuonPt1/F");
  all1Tree->Branch("dimuonPt1",&dimuonPt1,"dimuonPt1/F");
  matchTree->Branch("dimuonPt2",&dimuonPt2,"dimuonPt2/F");
  only2Tree->Branch("dimuonPt2",&dimuonPt2,"dimuonPt2/F");
  all2Tree->Branch("dimuonPt2",&dimuonPt2,"dimuonPt2/F");

  matchTree->Branch("muonLead_pt1",&muonLead_pt1,"muonLead_pt1/F");
  only1Tree->Branch("muonLead_pt1",&muonLead_pt1,"muonLead_pt1/F");
  all1Tree->Branch("muonLead_pt1",&muonLead_pt1,"muonLead_pt1/F");
  matchTree->Branch("muonLead_pt2",&muonLead_pt2,"muonLead_pt2/F");
  all2Tree->Branch("muonLead_pt2",&muonLead_pt2,"muonLead_pt2/F");

  matchTree->Branch("muonSub_pt1",&muonSub_pt1,"muonSub_pt1/F");
  only1Tree->Branch("muonSub_pt1",&muonSub_pt1,"muonSub_pt1/F");
  all1Tree->Branch("muonSub_pt1",&muonSub_pt1,"muonSub_pt1/F");
  matchTree->Branch("muonSub_pt2",&muonSub_pt2,"muonSub_pt2/F");
  only2Tree->Branch("muonSub_pt2",&muonSub_pt2,"muonSub_pt2/F");
  all2Tree->Branch("muonSub_pt2",&muonSub_pt2,"muonSub_pt2/F");

  matchTree->Branch("muonLead_eta1",&muonLead_eta1,"muonLead_eta1/F");
  only1Tree->Branch("muonLead_eta1",&muonLead_eta1,"muonLead_eta1/F");
  all1Tree->Branch("muonLead_eta1",&muonLead_eta1,"muonLead_eta1/F");
  matchTree->Branch("muonLead_eta2",&muonLead_eta2,"muonLead_eta2/F");
  only2Tree->Branch("muonLead_eta2",&muonLead_eta2,"muonLead_eta2/F");
  all2Tree->Branch("muonLead_eta2",&muonLead_eta2,"muonLead_eta2/F");

  matchTree->Branch("muonSub_eta1",&muonSub_eta1,"muonSub_eta1/F");
  only1Tree->Branch("muonSub_eta1",&muonSub_eta1,"muonSub_eta1/F");
  all1Tree->Branch("muonSub_eta1",&muonSub_eta1,"muonSub_eta1/F");
  matchTree->Branch("muonSub_eta2",&muonSub_eta2,"muonSub_eta2/F");
  only2Tree->Branch("muonSub_eta2",&muonSub_eta2,"muonSub_eta2/F");
  all2Tree->Branch("muonSub_eta2",&muonSub_eta2,"muonSub_eta2/F");

  matchTree->Branch("muonLead_phi1",&muonLead_phi1,"muonLead_phi1/F");
  only1Tree->Branch("muonLead_phi1",&muonLead_phi1,"muonLead_phi1/F");
  all1Tree->Branch("muonLead_phi1",&muonLead_phi1,"muonLead_phi1/F");
  matchTree->Branch("muonLead_phi2",&muonLead_phi2,"muonLead_phi2/F");
  only2Tree->Branch("muonLead_phi2",&muonLead_phi2,"muonLead_phi2/F");
  all2Tree->Branch("muonLead_phi2",&muonLead_phi2,"muonLead_phi2/F");

  matchTree->Branch("muonSub_phi1",&muonSub_phi1,"muonSub_phi1/F");
  only1Tree->Branch("muonSub_phi1",&muonSub_phi1,"muonSub_phi1/F");
  all1Tree->Branch("muonSub_phi1",&muonSub_phi1,"muonSub_phi1/F");
  matchTree->Branch("muonSub_phi2",&muonSub_phi2,"muonSub_phi2/F");
  only2Tree->Branch("muonSub_phi2",&muonSub_phi2,"muonSub_phi2/F");
  all2Tree->Branch("muonSub_phi2",&muonSub_phi2,"muonSub_phi2/F");

  matchTree->Branch("eventType1",&eventType1,"eventType1/I");
  only1Tree->Branch("eventType1",&eventType1,"eventType1/I");
  all1Tree->Branch("eventType1",&eventType1,"eventType1/I");
  matchTree->Branch("eventType2",&eventType2,"eventType2/I");
  only2Tree->Branch("eventType2",&eventType2,"eventType2/I");
  all2Tree->Branch("eventType2",&eventType2,"eventType2/I");

  matchTree->Branch("dimuonMass_noMuscle1",&dimuonMass_noMuscle1,"dimuonMass_noMuscle1/F");
  only1Tree->Branch("dimuonMass_noMuscle1",&dimuonMass_noMuscle1,"dimuonMass_noMuscle1/F");
  all1Tree->Branch("dimuonMass_noMuscle1",&dimuonMass_noMuscle1,"dimuonMass_noMuscle1/F");
  matchTree->Branch("dimuonMass_noMuscle2",&dimuonMass_noMuscle2,"dimuonMass_noMuscle2/F");
  only2Tree->Branch("dimuonMass_noMuscle2",&dimuonMass_noMuscle2,"dimuonMass_noMuscle2/F");
  all2Tree->Branch("dimuonMass_noMuscle2",&dimuonMass_noMuscle2,"dimuonMass_noMuscle2/F");

  matchTree->Branch("dimuonPt_noMuscle1",&dimuonPt_noMuscle1,"dimuonPt_noMuscle1/F");
  only1Tree->Branch("dimuonPt_noMuscle1",&dimuonPt_noMuscle1,"dimuonPt_noMuscle1/F");
  all1Tree->Branch("dimuonPt_noMuscle1",&dimuonPt_noMuscle1,"dimuonPt_noMuscle1/F");
  matchTree->Branch("dimuonPt_noMuscle2",&dimuonPt_noMuscle2,"dimuonPt_noMuscle2/F");
  only2Tree->Branch("dimuonPt_noMuscle2",&dimuonPt_noMuscle2,"dimuonPt_noMuscle2/F");
  all2Tree->Branch("dimuonPt_noMuscle2",&dimuonPt_noMuscle2,"dimuonPt_noMuscle2/F");

  matchTree->Branch("muonLead_pt_noMuscle1",&muonLead_pt_noMuscle1,"muonLead_pt_noMuscle1/F");
  only1Tree->Branch("muonLead_pt_noMuscle1",&muonLead_pt_noMuscle1,"muonLead_pt_noMuscle1/F");
  all1Tree->Branch("muonLead_pt_noMuscle1",&muonLead_pt_noMuscle1,"muonLead_pt_noMuscle1/F");
  matchTree->Branch("muonLead_pt_noMuscle2",&muonLead_pt_noMuscle2,"muonLead_pt_noMuscle2/F");
  all2Tree->Branch("muonLead_pt_noMuscle2",&muonLead_pt_noMuscle2,"muonLead_pt_noMuscle2/F");

  matchTree->Branch("muonSub_pt_noMuscle1",&muonSub_pt_noMuscle1,"muonSub_pt_noMuscle1/F");
  only1Tree->Branch("muonSub_pt_noMuscle1",&muonSub_pt_noMuscle1,"muonSub_pt_noMuscle1/F");
  all1Tree->Branch("muonSub_pt_noMuscle1",&muonSub_pt_noMuscle1,"muonSub_pt_noMuscle1/F");
  matchTree->Branch("muonSub_pt_noMuscle2",&muonSub_pt_noMuscle2,"muonSub_pt_noMuscle2/F");
  only2Tree->Branch("muonSub_pt_noMuscle2",&muonSub_pt_noMuscle2,"muonSub_pt_noMuscle2/F");
  all2Tree->Branch("muonSub_pt_noMuscle2",&muonSub_pt_noMuscle2,"muonSub_pt_noMuscle2/F");

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

  int nEvents1 = tree1->GetEntries();
  cout << "nEvents1: " << nEvents1 << endl;
  int nEvents2 = tree2->GetEntries();
  cout << "nEvents2: " << nEvents2 << endl;

//   result += " && muonLead_passPFRelIso && muonSub_passPFRelIso && ("#(muonLead_isHltMatched || muonSub_isHltMatched)"
//    result += " (muonLead_pt>25 && muonLead_isHltMatched) || "
//    result += " (muonSub_pt >25 &&  muonSub_isHltMatched) "
//    result += ")"


  for(int i=0; i<nEvents1;i++)
  {
    tree1->GetEntry(i);
    if(dimuonMass1<110.)
      continue;
    if(useTrkIso)
    {
      if (!(muonLead_passTrkRelIso1 && muonSub_passTrkRelIso1))
        continue;
    }
    else
    {
      if (!(muonLead_passPFRelIso1 && muonSub_passPFRelIso1))
        continue;
    }
    if(
        !(muonLead_pt1 > 25. && muonLead_isHltMatched1) 
        && !(muonSub_pt1 > 25 && muonSub_isHltMatched1)
    )
      continue;
    all1Tree->Fill();
    string eventStr = "";
    eventStr.appendAny(eventInfo_run1);
    eventStr.append(":");
    eventStr.appendAny(eventInfo_event1);
    if (events1.count(eventStr)==0)
    {
      events1[eventStr]=i;
    }
    else
    {
      cout << "Event "<<eventStr<<" found multiple times in first dataset!!\n";
    }
    runs1.insert(eventInfo_run1);
  }// end event1 loop

  for(int i=0; i<nEvents2;i++)
  {
    tree2->GetEntry(i);
    if(dimuonMass2<110.)
      continue;
    if(!useMuscleFitSecond) // default
    {
      if(useTrkIso)
      {
        if (!(muonLead_passTrkRelIso_noMuscle2 && muonSub_passTrkRelIso_noMuscle2))
          continue;
      }
      else
      {
        if (!(muonLead_passPFRelIso_noMuscle2 && muonSub_passPFRelIso_noMuscle2))
          continue;
      }
      if(
          !(muonLead_pt_noMuscle2 > 25. && muonLead_isHltMatched2) 
          && !(muonSub_pt_noMuscle2 > 25 && muonSub_isHltMatched2)
      )
        continue;
    }
    else
    {
      if(useTrkIso)
      {
        if (!(muonLead_passTrkRelIso2 && muonSub_passTrkRelIso2))
          continue;
      }
      else
      {
        if (!(muonLead_passPFRelIso2 && muonSub_passPFRelIso2))
          continue;
      }
      if(
          !(muonLead_pt2 > 25. && muonLead_isHltMatched2) 
          && !(muonSub_pt2 > 25 && muonSub_isHltMatched2)
      )
        continue;
    }
    all2Tree->Fill();
    string eventStr = "";
    eventStr.appendAny(eventInfo_run2);
    eventStr.append(":");
    eventStr.appendAny(eventInfo_event2);
    if (events2.count(eventStr)==0)
    {
      events2[eventStr]=i;
    }
    else
    {
      cout << "Event "<<eventStr<<" found multiple times in second dataset!!\n";
    }
    runs2.insert(eventInfo_run2);
  }// end event2 loop

  unsigned nevents1=0;
  for(events1Itr=events1.begin();events1Itr!=events1.end();events1Itr++)
  {
    nevents1++;
    if(events2.count(events1Itr->first)>0)
    {
      tree1->GetEntry(events1Itr->second);
      tree2->GetEntry(events2[events1Itr->first]);
      if(eventInfo_event1 != eventInfo_event2 || eventInfo_run1 != eventInfo_run2)
      {
        cout << "Error: Events don't match!!\n";
        continue;
      }
      matchTree->Fill();
    }
    else
    {
      tree1->GetEntry(events1Itr->second);
      only1Tree->Fill();
    }
  }

  unsigned nevents2=0;
  for(events2Itr=events2.begin();events2Itr!=events2.end();events2Itr++)
  {
    nevents2++;
    if(events1.count(events2Itr->first)==0)
    {
      tree2->GetEntry(events2Itr->second);
      only2Tree->Fill();
    }
  }

  matchTree->Write();
  only1Tree->Write();
  only2Tree->Write();
  all1Tree->Write();
  all2Tree->Write();

  cout << "\nMatches: " << matchTree->GetEntries()<<"\n";
  cout << "Only 1: " << only1Tree->GetEntries()<<"\n";
  cout << "Only 2: " << only2Tree->GetEntries()<<"\n";

  cout << "\nevents1.size(): "<< events1.size() << endl;
  cout << "events2.size(): "<< events2.size() << endl;

  /*
  events1Itr=events1.begin();
  events2Itr=events2.begin();
  while(true)
  {
    string str1;
    string str2;
    bool end1=false;
    bool end2=false;
    if (events1Itr!=events1.end())
    {
      str1 = events1Itr->first;
      events1Itr++;
    }
    else
      end1=true; 
    if (events2Itr!=events2.end())
    {
      str2 = events2Itr->first;
      events2Itr++;
    }
    else
      end2=true; 
    
    if(end1 && end2)
      break;

    cout << setw(20)<< str1 << setw(20)<< str2<<endl;
  }

  runs1Itr=runs1.begin();
  runs2Itr=runs2.begin();
  while(true)
  {
    string str1;
    string str2;
    bool end1=false;
    bool end2=false;
    if (runs1Itr!=runs1.end())
    {
      str1.appendAny(*runs1Itr);
      runs1Itr++;
    }
    else
      end1=true; 
    if (runs2Itr!=runs2.end())
    {
      str2.appendAny(*runs2Itr);
      runs2Itr++;
    }
    else
      end2=true; 
    
    if(end1 && end2)
      break;

    cout << setw(20)<< str1 << setw(20)<< str2<<endl;
  }
  */

  outfile->Write();
  outfile->Close();

  return 0;
}
