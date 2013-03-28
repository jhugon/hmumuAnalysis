#include "MEKD_Wrapper.h"
#include <iostream>

MEKD_Wrapper::MEKD_Wrapper(double collisionEnergy /*in TeV*/ )
  : _collisionEnergy(collisionEnergy), _mekd()
{
  _mekd.Sqrt_s = collisionEnergy*1000.; //in GeV
  _mekd.Use_PDF_w_pT0 = false;
  _mekd.Resonance_decay_mode = "2mu";
  _signalProcessName = "ggSpin0SMH";
  _backgroundProcessName = "DY";

//  _mekd.Debug_Mode = true;
//  _mekd.Warning_Mode = true;

//  std::cout << "MEKD Original Card and PDF:\n";
//  std::cout << _mekd.Parameter_file << std::endl;
//  std::cout << _mekd.PDF_file << std::endl;

  _mekd.Parameter_file = "mekd/src/Cards/param_card.dat";
  _mekd.PDF_file = "mekd/src/PDFTables/cteq6l.pdt";

//  std::cout << "MEKD Changed Card and PDF:\n";
//  std::cout << _mekd.Parameter_file << std::endl;
//  std::cout << _mekd.PDF_file << std::endl;
}

int
MEKD_Wrapper::getKD(const TLorentzVector& mu1, 
                    const TLorentzVector& mu2, 
                    int mu1Charge,
                    float& kd, float& sigME, float& bakME)
{
  using namespace std;
  kd = -1.;
  sigME = -1.;
  bakME = -1.;
  _mekd.p3 = 0;
  _mekd.id3 = 0;
  _mekd.p4 = 0;
  _mekd.id4 = 0;
  _mekd.p5 = 0;
  _mekd.id5 = 0;

  if (mu1Charge < 0.0)
  {
    _mekd.id1=13;
    _mekd.id2=-13;
  }
  else
  {
    _mekd.id1=-13;
    _mekd.id2=13;
  }

  double p1[4];
  double p2[4];

  p1[0] = mu1.E();
  p1[1] = mu1.Px();
  p1[2] = mu1.Py();
  p1[3] = mu1.Pz();
  _mekd.p1 = p1;

  p2[0] = mu2.E();
  p2[1] = mu2.Px();
  p2[2] = mu2.Py();
  p2[3] = mu2.Pz();
  _mekd.p2 = p2;
  
  int result = _mekd.Run_MEKD_MG(_signalProcessName);
  sigME = _mekd.Signal_ME;
  result += 10*_mekd.Run_MEKD_MG(_backgroundProcessName);
  bakME = _mekd.Signal_ME;

  kd = log(sigME/bakME);

  return result;
}
