#ifndef helpers_h
#define helpers_h
#include "DataFormats.h"
#include <algorithm>
#include "TRandom.h"
#include "TLorentzVector.h"

enum PUJetID
{
  puJetTight  = 0,
  puJetMedium = 1,
  puJetLoose  = 2
};

enum SelectionCodes
{
  notSelected_Code = 0,

  vbfPresel_Code = 1,
  vbfBDTCut_Code = 2,
  incPresel_Code = 4,
  incBDTCut_Code = 8,
  isBB_Code = 16,
  isBO_Code = 32,
  isBE_Code = 64,
  isOO_Code = 128,
  isOE_Code = 256,
  isEE_Code = 512,
  isNotBB_Code = 1024,
};

float getPFRelIso(_MuonInfo& muon);
float getTrkRelIso(_MuonInfo& muon);

bool isKinTight_2012(_MuonInfo& muon);
bool isKinTight_2012_noIso(_MuonInfo& muon);

bool isKinTight_2011(_MuonInfo& muon);
bool isKinTight_2011_noIso(_MuonInfo& muon);

bool passPUJetID(int flag, PUJetID desiredLevel);

float smearMC(float trueVal, float recoVal, float calib, float smearRatio,TRandom random, bool debug=false);

bool isHltMatched(_MuonInfo& muon1, _MuonInfo& muon2, std::vector<int> allowedPaths);
bool isHltMatched(_MuonInfo& muon1, std::vector<int> allowedPaths);

float resolutionBias(float eta);
float jerCorr   (float ptold, float oldgenpt, float etaold);
float corrPtUp  (float ptold, float oldgenpt, float etaold);
float corrPtDown(float ptold, float oldgenpt, float etaold);

int whichSelection(_MuonInfo& mu1, _MuonInfo& mu2,
                   std::vector<int> allowedHLTPaths,
                   std::string& runPeriod,
                   _PFJetInfo jets,
                   bool passIncBDTCut,
                   bool passVBFBDTCut,
                   double sigmasJEC=0,
                   double jerUncertainty=0);

int whichSelection(_MuonInfo& mu1, _MuonInfo& mu2,
                   std::string& runPeriod,
                   _PFJetInfo jets,
                   bool passIncBDTCut,
                   bool passVBFBDTCut);

void setStyle();

#endif
