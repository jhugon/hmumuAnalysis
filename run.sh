#!/bin/bash

nice scons -j4

nice ./analyzer ggHmumu125.root /data/uftrig01b/jhugon/hmumu/dataPU/ggHmumu125.root
nice ./analyzer vbfHmumu125.root /data/uftrig01b/jhugon/hmumu/dataPU/vbfHmumu125.root

#nice ./analyzer DYJetsToLL.root /data/uftrig01b/digiovan/root/CMSSW_5_2_5/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/*.root
