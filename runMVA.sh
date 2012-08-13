#!/bin/bash

nice scons

nice ./trainingTreeMaker signalTree.root ../dataPU/ggHmumu125.root
nice ./trainingTreeMaker backgroundTree.root /data/uftrig01b/digiovan/root/CMSSW_5_2_5/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/*.root

nice root -b -q -x TMVAClassification.C++
