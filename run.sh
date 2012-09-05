#!/bin/bash

nice scons -j4

nice ./analyzer /data/uftrig01b/jhugon/hmumu/dataPU/ggHmumu125.root  ggHmumu125.root --trainingTree signalTreeGG.root
nice ./analyzer /data/uftrig01b/jhugon/hmumu/dataPU/vbfHmumu125.root vbfHmumu125.root --trainingTree signalTreeVBF.root
nice ./analyzer /data/uftrig01b/jhugon/hmumu/backgroundNtuples/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root  DYJetsToLL.root --trainingTree backgroundTreeDY.root
nice ./analyzer /data/uftrig01b/jhugon/hmumu/backgroundNtuples/NtuplesMCTTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1.root ttbar.root --trainingTree backgroundTreeTT.root

./mvaTrain inclusive.cfg >& logMVA

# Run with full MVA
nice ./analyzer /data/uftrig01b/jhugon/hmumu/dataPU/ggHmumu125.root  ggHmumu125.root
nice ./analyzer /data/uftrig01b/jhugon/hmumu/dataPU/vbfHmumu125.root vbfHmumu125.root
nice ./analyzer /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch1/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root  DYJetsToLL.root

nice ./analyzer /data/uftrig01b/digiovan/root/higgs/CMSSW_5_2_6/NtuplesMCTTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/*.root ttbar.root
