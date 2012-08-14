#!/bin/bash

nice scons

nice ./trainingTreeMaker signalTreeInc.root ../dataPU/ggHmumu125.root
nice ./trainingTreeMaker signalTreeVBF.root ../dataPU/vbfHmumu125.root
nice ./trainingTreeMaker backgroundTree.root /data/uftrig01b/digiovan/root/CMSSW_5_2_5/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/*.root

rm -rf weights

nice root -b -q -x TMVAClassificationVBF.C++ >& logMVAVBF
rm -rf weightsVBF/
mv weights weightsVBF
nice root -b -q -x TMVAClassificationInc.C++ >& logMVAInc
rm -rf weightsMuonOnly
mv weights weightsMuonOnly
