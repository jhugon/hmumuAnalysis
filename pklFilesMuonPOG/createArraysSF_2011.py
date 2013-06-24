#!/usr/bin/env python

import pickle

# open the files
file = open('MuonEfficiencies2011_44X.pkl', 'r')

# dictionary
dict = pickle.load( file )

# defines the eta/pt ranges
etaRange2011 = ['-2.1_-1.6',
            '-1.6_-1.2',
            '-1.2_-0.9',
            '-0.9_-0.6',
            '-0.6_-0.3',
            '-0.3_-0.2',
            '-0.2_0.2',
            '0.2_0.3',
            '0.3_0.6',
            '0.6_0.9',
            '0.9_1.2',                  
            '1.2_1.6',
            '1.6_2.1']

ptRange2011 = ['10_20',
           '20_30',
           '30_40',
           '40_50',
           '50_60',
           '60_80',
           '80_250']


print "const int etaBins2011 = 13;"
print "const int ptBins2011 = 7;"
print " "
print "double etaRange2011[etaBins2011+1] = {-2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2,"
print "                              0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1};"
print " "
print "double ptRange2011[ptBins2011+1] = {10, 20, 30, 40, 50, 60, 80, 250};"
print " "
print " "

print "// 2011A ID/Iso" 
print "double MuonIdTight_2011A[etaBins2011] = {"
for eta in etaRange2011:
    print dict['TIGHT_nL8_2011A']['eta_pt>20'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta12_2011A[ptBins2011] = {"
for pt in ptRange2011:
    print dict['combRelPFISO12_2011A']['pt_abseta<1.2'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta12to21_2011A[ptBins2011] = {"
for pt in ptRange2011:
    print dict['combRelPFISO12_2011A']['pt_abseta>1.2'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"





print "// 2011B ID/Iso" 
print "double MuonIdTight_2011B[etaBins2011] = {"
for eta in etaRange2011:
    print dict['TIGHT_nL8_2011B']['eta_pt>20'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta12_2011B[ptBins2011] = {"
for pt in ptRange2011:
    print dict['combRelPFISO12_2011B']['pt_abseta<1.2'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta12to21_2011B[ptBins2011] = {"
for pt in ptRange2011:
    print dict['combRelPFISO12_2011B']['pt_abseta>1.2'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"

