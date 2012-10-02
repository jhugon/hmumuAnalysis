#ifndef helpers_h
#define helpers_h
#include "DataFormats.h"
#include <algorithm>
#include "TRandom.h"

enum PUJetID
{
  puJetTight  = 0,
  puJetMedium = 1,
  puJetLoose  = 2
};

float getRelIso(_MuonInfo& muon);

bool isKinTight_2012(_MuonInfo& muon);

bool passPUJetID(int flag, PUJetID desiredLevel);

float smearMC(float trueVal, float recoVal, float calib, float smearRatio,TRandom random, bool debug=false);

#endif
