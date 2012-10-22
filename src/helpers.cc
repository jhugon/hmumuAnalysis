#include "helpers.h"
#include <cmath>
#include <iostream>

float getRelIso(_MuonInfo& muon)
{
  // tracker iso cut
  //if (muon.trackIsoSumPt/muon.pt>0.1) return isKinTight_2012;  

  float result = muon.sumChargedHadronPtR04 + 
    std::max(0.0,muon.sumNeutralHadronEtR04 + muon.sumPhotonEtR04 - 0.5*muon.sumPUPtR04);

  return result/muon.pt;
}

bool isKinTight_2012(_MuonInfo& muon) 
{

  bool isKinTight_2012=false;

  if (!muon.isGlobal)  return isKinTight_2012;
  if (!muon.isPFMuon) return isKinTight_2012;

  // acceptance cuts
  if (muon.pt < 25)         return isKinTight_2012; // pt cut
  if (fabs(muon.eta) > 2.1) return isKinTight_2012; // eta cut

  // kinematic cuts
  if (muon.numTrackerLayers < 6) return isKinTight_2012; // # hits in tracker

  if(getRelIso(muon) > 0.12)
      return isKinTight_2012;

  if (fabs(muon.d0_PV) > 0.2) return isKinTight_2012;
  if (fabs(muon.dz_PV) > 0.5) return isKinTight_2012;

  if ( muon.numValidMuonHits  < 1 ) return isKinTight_2012;
  if ( muon.numValidPixelHits < 1 ) return isKinTight_2012;
  if ( muon.numOfMatchedStations < 2 ) return isKinTight_2012;
  if ( muon.normChiSquare > 10)     return isKinTight_2012;

  isKinTight_2012=true;
  return isKinTight_2012;
}

bool passPUJetID(int flag, PUJetID desiredLevel)
{
 return ( flag & (1 << desiredLevel) ) != 0;
}

float smearMC(float trueVal, float recoVal, float calib, float smearRatio,TRandom random, bool debug)
{
  if (trueVal > 0.0)
  {
    float uncalibratedReco = recoVal+calib;
    float dif = uncalibratedReco-trueVal;
    dif *= smearRatio;
    float result = trueVal+dif;
    if (debug)
      std::cout << "true: " << trueVal << " \t\tuncalReco: " << uncalibratedReco << " \t\trecoCand: " << recoVal << " \t\tfinal: " << result << std::endl;
    return result;
  }
  else
  {
    return recoVal+calib;
  }
}

bool isHltMatched(_MuonInfo& muon1, _MuonInfo& muon2, std::vector<int> allowedPaths)
{
  std::vector<int>::const_iterator path = allowedPaths.begin();
  std::vector<int>::const_iterator endPath = allowedPaths.end();
  for(;path != endPath;path++)
  {
    //std::cout << "muon1 index "<<*path<<": "<<muon1.isHltMatched[*path]<<std::endl;
    //std::cout << "muon2 index "<<*path<<": "<<muon2.isHltMatched[*path]<<std::endl;
    if(muon1.isHltMatched[*path]==1)
        return true;
    if(muon2.isHltMatched[*path]==1)
        return true;
  }
  return false;
}

bool isHltMatched(_MuonInfo& muon1, std::vector<int> allowedPaths)
{
  std::vector<int>::const_iterator path = allowedPaths.begin();
  std::vector<int>::const_iterator endPath = allowedPaths.end();
  for(;path != endPath;path++)
  {
    //std::cout << "muon1 index "<<*path<<": "<<muon1.isHltMatched[*path]<<std::endl;
    if(muon1.isHltMatched[*path]==1)
        return true;
  }
  return false;
}
