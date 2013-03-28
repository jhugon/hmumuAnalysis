//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph 5 v. 1.5.7, 2013-01-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef HelAmps_HEF_MEKD2_spin2_2l_H
#define HelAmps_HEF_MEKD2_spin2_2l_H

#include <cmath> 
#include <complex> 

using namespace std; 

namespace MG5_HEF_MEKD2_spin2_2l
{
double Sgn(double e, double f); 

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void VVT5_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT13_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFV2_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[]);

void FFV2_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[]);

void FFT3_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFT4_0(complex<double> F1[], complex<double> F2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void FFT4_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFT2_0(complex<double> F1[], complex<double> F2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT10_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT12_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT9_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT4_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFT1_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFT1_2_3_4_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, complex<double> COUP4,
    double M3, double W3, complex<double> T3[]);

void VVT11_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT1_10_11_12_13_3_5_7_8_9_3(complex<double> V1[], complex<double> V2[],
    complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> COUP5, complex<double> COUP6,
    complex<double> COUP7, complex<double> COUP8, complex<double> COUP9,
    complex<double> COUP10, double M3, double W3, complex<double> T3[]);

void FFT3_0(complex<double> F1[], complex<double> F2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void VVT7_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT8_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void VVT3_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFT2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> T3[]);

void FFT1_0(complex<double> F1[], complex<double> F2[], complex<double> T3[],
    complex<double> COUP, complex<double> & vertex);

void FFT1_2_3_4_0(complex<double> F1[], complex<double> F2[], complex<double>
    T3[], complex<double> COUP1, complex<double> COUP2, complex<double> COUP3,
    complex<double> COUP4, complex<double> & vertex);

}  // end namespace MG5_HEF_MEKD

#endif  // HelAmps_HEF_MEKD2_H
