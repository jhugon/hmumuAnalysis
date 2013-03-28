#ifndef MEKD_WRAPPER_h
#define MEKD_WRAPPER_h

#include "MEKD_MG.h"
#include "TLorentzVector.h"

class MEKD_Wrapper
{
  public:
    MEKD_Wrapper(double collisionEnergy /*in TeV*/);
    int getKD(const TLorentzVector& mu1, const TLorentzVector& mu2, 
                            int mu1Charge, 
                            float& kd, float& sigME, float& bakME);

  private:
    double _collisionEnergy;
    MEKD_MG _mekd;
    std::string _signalProcessName;
    std::string _backgroundProcessName;
};

#endif
