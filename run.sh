#!/bin/bash

nice scons -j4

#nice ./analyzer ggHmumu.root /data/uftrig01b/jhugon/hmumu/dataNoPUBadNtuples/ggHmumu125-NoPU.root
#nice ./analyzer vbfHmumu.root /data/uftrig01b/jhugon/hmumu/dataNoPUBadNtuples/vbfHmumu125-NoPU.root

nice ./analyzer ggHmumu125.root /data/uftrig01b/jhugon/hmumu/dataPU/ggHmumu125.root
nice ./analyzer vbfHmumu150.root /data/uftrig01b/jhugon/hmumu/dataPU/vbfHmumu150.root

nice ./analyzer DYJetsToLL.root /data/uftrig01b/digiovan/root/CMSSW_5_2_5/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2/*.root
