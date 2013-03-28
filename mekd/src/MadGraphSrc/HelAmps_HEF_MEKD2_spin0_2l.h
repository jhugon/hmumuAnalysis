//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph 5 v. 1.5.7, 2013-01-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef HelAmps_HEF_MEKD2_spin0_2l_H
#define HelAmps_HEF_MEKD2_spin0_2l_H

#include <cmath> 
#include <complex> 

using namespace std; 

namespace MG5_HEF_MEKD2_spin0_2l
{
double Sgn(double e, double f); 

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void FFS2_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);

void VVS3_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);
void VVS3_4_5_3(complex<double> V1[], complex<double> V2[], complex<double>
    COUP1, complex<double> COUP2, complex<double> COUP3, double M3, double W3,
    complex<double> S3[]);

void FFS1_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);
void FFS1_2_3(complex<double> F1[], complex<double> F2[], complex<double>
    COUP1, complex<double> COUP2, double M3, double W3, complex<double> S3[]);

void VVS5_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);

void VVS1_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);

void VVS4_3(complex<double> V1[], complex<double> V2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);

void FFS1_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);
void FFS1_2_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex);

void FFS2_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> S3[]);

}  // end namespace MG5_HEF_MEKD

#endif  // HelAmps_HEF_MEKD2_H
