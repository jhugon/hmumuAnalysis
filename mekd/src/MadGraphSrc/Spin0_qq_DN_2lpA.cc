//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.7, 2013-01-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "Spin0_qq_DN_2lpA.h"
#include "HelAmps_HEF_MEKD2_spin0_2lpA.h"

using namespace MG5_HEF_MEKD2_spin0_2lpA; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: s s~ > h > mu- mu+ a / zp xg WEIGHTED=4

//--------------------------------------------------------------------------
// Initialize process.

void Spin0_qq_DN_2lpA::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_HEF_MEKD::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings();
  // Set external particle masses for this matrix element
  mME.push_back(pars->MS); 
  mME.push_back(pars->MS); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->MM); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[1]; 
}

//--------------------------------------------------------------------------
// Update process.

void Spin0_qq_DN_2lpA::updateProc(SLHAReader slha) 
{
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  
  // Set external particle masses for this matrix element
  mME[0]=(pars->MS);
  mME[1]=(pars->MS);
  mME[2]=(pars->MM);
  mME[3]=(pars->MM);
  mME[4]=(pars->ZERO);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void Spin0_qq_DN_2lpA::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings();

  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 32; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
//   std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1}, {-1,
      -1, -1, -1, 1}, {-1, -1, -1, 1, -1}, {-1, -1, -1, 1, 1}, {-1, -1, 1, -1,
      -1}, {-1, -1, 1, -1, 1}, {-1, -1, 1, 1, -1}, {-1, -1, 1, 1, 1}, {-1, 1,
      -1, -1, -1}, {-1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1},
      {-1, 1, 1, -1, -1}, {-1, 1, 1, -1, 1}, {-1, 1, 1, 1, -1}, {-1, 1, 1, 1,
      1}, {1, -1, -1, -1, -1}, {1, -1, -1, -1, 1}, {1, -1, -1, 1, -1}, {1, -1,
      -1, 1, 1}, {1, -1, 1, -1, -1}, {1, -1, 1, -1, 1}, {1, -1, 1, 1, -1}, {1,
      -1, 1, 1, 1}, {1, 1, -1, -1, -1}, {1, 1, -1, -1, 1}, {1, 1, -1, 1, -1},
      {1, 1, -1, 1, 1}, {1, 1, 1, -1, -1}, {1, 1, 1, -1, 1}, {1, 1, 1, 1, -1},
      {1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {36, 36}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_ssx_h_mummupa_no_zpxg(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_ssx_h_mummupa_no_zpxg(); 
        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_ssx_h_mummupa_no_zpxg(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_ssx_h_mummupa_no_zpxg(); 
      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double Spin0_qq_DN_2lpA::sigmaHat() 
{
  // Select between the different processes
  if(id1 == -3 && id2 == 3)
  {
    // Add matrix elements for processes with beams (-3, 3)
    return matrix_element[1]; 
  }
  else if(id1 == 3 && id2 == -3)
  {
    // Add matrix elements for processes with beams (3, -3)
    return matrix_element[0]; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void Spin0_qq_DN_2lpA::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//   int i, j; 

  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  FFS1_2_3(w[0], w[1], pars->HEF_MEKD2_GC_137, pars->HEF_MEKD2_GC_136, pars->MH, pars->WH, w[5]); 
  FFV2_1(w[2], w[4], pars->HEF_MEKD2_GC_5, pars->MM, pars->ZERO, w[6]); 
  FFV2_2(w[3], w[4], pars->HEF_MEKD2_GC_5, pars->MM, pars->ZERO, w[7]); 
  FFV2_2(w[0], w[4], pars->HEF_MEKD2_GC_3, pars->MS, pars->ZERO, w[8]); 
  FFS1_2_3(w[3], w[2], pars->HEF_MEKD2_GC_127, pars->HEF_MEKD2_GC_126, pars->MH, pars->WH, w[9]); 
  FFV2_1(w[1], w[4], pars->HEF_MEKD2_GC_3, pars->MS, pars->ZERO, w[10]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFS1_2_0(w[3], w[6], w[5], pars->HEF_MEKD2_GC_127, pars->HEF_MEKD2_GC_126, amp[0]); 
  FFS1_2_0(w[7], w[2], w[5], pars->HEF_MEKD2_GC_127, pars->HEF_MEKD2_GC_126, amp[1]); 
  FFS1_2_0(w[8], w[1], w[9], pars->HEF_MEKD2_GC_137, pars->HEF_MEKD2_GC_136, amp[2]); 
  FFS1_2_0(w[0], w[10], w[9], pars->HEF_MEKD2_GC_137, pars->HEF_MEKD2_GC_136, amp[3]); 

}
double Spin0_qq_DN_2lpA::matrix_ssx_h_mummupa_no_zpxg() 
{
  int i, j; 
  // Local variables
//   const int ngraphs = 4; 
  const int ncolor = 1; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{3}}; 

  // Calculate color flows
  jamp[0] = -amp[0] - amp[1] - amp[2] - amp[3]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



