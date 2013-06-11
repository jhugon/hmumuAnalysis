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

  vector<string> outFileNames;
  vector<string>::const_iterator outFileItr;
  if (optionMap.count("outfilename")>0)
  {
     outFileNames = optionMap["outfilename"].as<vector<string> >();

  }
  else
  {
     cout << "Error: outfilename  required, exiting." << endl;
     return 1;
  }

  for(outFileItr=outFileNames.begin();outFileItr!=outFileNames.end();outFileItr++)
  {
    TChain* tree1 = new TChain("outtree");
    tree1->AddFile(outFileItr->c_str());
  
    set<int> runs;
    set<int>::const_iterator runsItr;
  
    int eventInfo_run;
    tree1->SetBranchAddress("eventInfo_run", &eventInfo_run);
  
  
    int nEvents1 = tree1->GetEntries();
    //cout << "nEvents1: " << nEvents1 << endl;

    cout << *outFileItr;
  
    for(int i=0; i<nEvents1;i++)
    {
      tree1->GetEntry(i);
      runs.insert(eventInfo_run);
    }// end event1 loop

    delete tree1;

    for(runsItr=runs.begin();runsItr!=runs.end();runsItr++)
      cout << " "<< (*runsItr);
    cout << endl;
  }

  return 0;
}
