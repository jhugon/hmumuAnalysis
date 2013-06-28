#!/usr/bin/env python

import pickle

# open the files
# this is for ID and isolation
f_2012ABCD_ID_ISO = open('Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.pkl', 'r')

# these are old files for ID and isolation, but they are needed
# for the TRIGGER information!
f_2012AB = open('MuonEfficiencies_Run_2012A_2012_B_53X.pkl', 'r')
f_2012C  = open('MuonEfficiencies_Run_2012C_53X.pkl', 'r')
f_2012D_Trigger = open('TriggerMuonEfficiencies_Run_2012D_53X.pkl', 'r')  

# dictionary
dict_2012ABCD_ID_ISO = pickle.load( f_2012ABCD_ID_ISO )

dict_2012AB          = pickle.load( f_2012AB )
dict_2012C           = pickle.load( f_2012C  )
dict_2012D_Trigger   = pickle.load( f_2012D_Trigger )



# defines the eta/pt ranges
etaRange = ['-2.1_-1.6',
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

ptRange = ['10_20',
           '20_25',
           '25_30',
           '30_35',
           '35_40',
           '40_50',
           '50_60',
           '60_90',
           '90_140',
           '140_300',
           '300_500']


print "const int etaBins = 13;"
print "const int ptBins = 11;"
print " "
print "double etaRange[etaBins+1] = {-2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2,"
print "                              0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1};"
print " "
print "double ptRange[ptBins+1] = {10, 20, 25, 30, 35, 40, 50, 60, 90, 140, 300, 500};"
print " "
print " "

print "// 2012ABCD ID/Iso" 
print "double MuonIdTight_2012ABCD[etaBins] = {"
for eta in etaRange:
    print dict_2012ABCD_ID_ISO['Tight']['etapt20-500_2012ABCD'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"

print "double PFIsoIddB_eta09_2012ABCD[ptBins] = {"
for pt in ptRange:             
    print dict_2012ABCD_ID_ISO['combRelIsoPF04dBeta<012_Tight']['ptabseta<0.9_2012ABCD'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"

print "double PFIsoIddB_eta09to12_2012ABCD[ptBins] = {"
for pt in ptRange:
    print dict_2012ABCD_ID_ISO['combRelIsoPF04dBeta<012_Tight']['ptabseta0.9-1.2_2012ABCD'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta12to21_2012ABCD[ptBins] = {"
for pt in ptRange:
    print dict_2012ABCD_ID_ISO['combRelIsoPF04dBeta<012_Tight']['ptabseta1.2-2.1_2012ABCD'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"

##
## T R I G G E R
##
print "// Trigger efficiencies split by run periods"

print "// 2012 A "
print "double TrigId_2012A[etaBins] = {"
for eta in etaRange:
    print dict_2012AB['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012A'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"

print "double effData_TrigId_2012A[etaBins] = {"
for eta in etaRange:
    print dict_2012AB['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012A'][eta]['data']['efficiency'], ", "
print "};\n\n"

print "double effMC_TrigId_2012A[etaBins] = {"
for eta in etaRange:
    print dict_2012AB['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012A'][eta]['mc']['efficiency'], ", "
print "};\n\n"




print "// 2012 B "
print "double TrigId_2012B[etaBins] = {"
for eta in etaRange:
    print dict_2012AB['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012B'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"

print "double effData_TrigId_2012B[etaBins] = {"
for eta in etaRange:
    print dict_2012AB['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012B'][eta]['data']['efficiency'], ", "
print "};\n\n"

print "double effMC_TrigId_2012B[etaBins] = {"
for eta in etaRange:
    print dict_2012AB['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012B'][eta]['mc']['efficiency'], ", "
print "};\n\n"




print "// 2012 C "
print "double TrigId_2012C[etaBins] = {"
for eta in etaRange:
    print dict_2012C['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"

print "double effData_TrigId_2012C[etaBins] = {"
for eta in etaRange:
    print dict_2012C['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500'][eta]['data']['efficiency'], ", "
print "};\n\n"

print "double effMC_TrigId_2012C[etaBins] = {"
for eta in etaRange:
    print dict_2012C['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500'][eta]['mc']['efficiency'], ", "
print "};\n\n"




print "// 2012 D "
print "double TrigId_2012D[etaBins] = {"
for eta in etaRange:
    print dict_2012D_Trigger['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012D'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double effData_TrigId_2012D[etaBins] = {"
for eta in etaRange:
    print dict_2012D_Trigger['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012D'][eta]['data']['efficiency'], ", "
print "};\n\n"

print "double effMC_TrigId_2012D[etaBins] = {"
for eta in etaRange:
    print dict_2012D_Trigger['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012D'][eta]['mc']['efficiency'], ", "
print "};\n\n"



