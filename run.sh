#!/bin/bash

scons

#./analyzer ggHmumu.root /data/uftrig01b/jhugon/hmumu/dataNoPU/ggHmumu125-NoPU.root
#./analyzer vbfHmumu.root /data/uftrig01b/jhugon/hmumu/dataNoPU/vbfHmumu125-NoPU.root
#./analyzer DYJetsToLL.root /data/uftrig01b/digiovan/root/CMSSW_5_2_5/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/*.root

./analyzer vbfHmumu.root /data/uftrig01b/jhugon/hmumu/dataPU/vbfHmumu150.root
