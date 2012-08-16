#!/bin/bash

nice scons

#nice ./trainingTreeMaker signalTreeInc.root ../dataPU/ggHmumu125.root
#nice ./trainingTreeMaker signalTreeVBF.root ../dataPU/vbfHmumu125.root
#nice ./trainingTreeMaker backgroundTreeDY.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch1/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root
#nice ./trainingTreeMaker backgroundTreeTT.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_2_6/NtuplesMCTTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/*.root

rm -rf weights
nice root -b -q -x TMVAClassificationVBF.C++ >& logMVAVBF
rm -rf weightsVBF/
mv weights weightsVBF
nice root -b -q -x TMVAClassificationInc.C++ >& logMVAInc
rm -rf weightsMuonOnly
mv weights weightsMuonOnly
