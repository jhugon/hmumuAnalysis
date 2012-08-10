#ifndef helpers_h
#define helpers_h
#include "DataFormats.h"
#include <algorithm>

bool isKinTight_2012(_MuonInfo& muon) 
{

  bool isKinTight_2012=false;

  if (!muon.isGlobal)  return isKinTight_2012;
  if (!muon.isPFMuon) return isKinTight_2012;

  // acceptance cuts
  if (muon.pt < 20)         return isKinTight_2012; // pt cut
  if (fabs(muon.eta) > 2.1) return isKinTight_2012; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return isKinTight_2012; // # hits in tracker

  // tracker iso cut
  //if (muon.trackIsoSumPt/muon.pt>0.1) return isKinTight_2012;  

  if ( ( muon.sumChargedHadronPtR04 + 
    std::max(0.0,muon.sumNeutralHadronEtR04 + muon.sumPhotonEtR04 - 0.5*muon.sumPUPtR04) ) / muon.pt > 0.12) 
      return isKinTight_2012;

  // beam spot cut
  if (fabs(muon.d0) > 0.2) return isKinTight_2012;
  // Will need dz cut when we have ntuplizer correct!!!

  if ( muon.numValidMuonHits  < 1 ) return isKinTight_2012;
  if ( muon.numValidPixelHits < 1 ) return isKinTight_2012;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_2012;
  if ( muon.normChiSquare > 10)     return isKinTight_2012;

  isKinTight_2012=true;
  return isKinTight_2012;
}

enum PUJetID
{
  kTight  = 0,
  kMedium = 1,
  kLoose  = 2
};

bool passPUJetID(int flag, puJetId desiredLevel)
{
 return ( flag & (1 << desiredLevel) ) != 0;
}

#endif
