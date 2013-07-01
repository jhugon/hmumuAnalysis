#!/bin/bash

nice scons -j4

#OPTIONS=" -m 1000"

##### Standard Samples

echo "Running 8 TeV Samples..."
DIR=/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/

nice ./analyzer ggHmumu115_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-115_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-115.root -r 8TeV $OPTIONS >& log_GluGlu-115.log & 
nice ./analyzer ggHmumu125_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-125_forXcheck.root -r 8TeV $OPTIONS >& log_GluGlu-125.log &
nice ./analyzer ggHmumu135_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-135_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-135.root -r 8TeV $OPTIONS >& log_GluGlu-135.log &
nice ./analyzer ggHmumu150_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-150_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-150.root -r 8TeV $OPTIONS >& log_GluGlu-150.log 

nice ./analyzer vbfHmumu115_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-115_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-115.root -r 8TeV $OPTIONS >& log_VBF-115.log & 
nice ./analyzer vbfHmumu125_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-125_forXcheck.root -r 8TeV $OPTIONS >& log_VBF-125.log &
nice ./analyzer vbfHmumu135_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-135_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-135.root -r 8TeV $OPTIONS >& log_VBF-135.log &
nice ./analyzer vbfHmumu150_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-150_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-150.root -r 8TeV $OPTIONS >& log_VBF-150.log &

wait

echo "Running 7 TeV Samples..."
DIR=/data/uftrig01b/digiovan/root/higgs/CMSSW_4_4_5/V00-01-10/

nice ./analyzer ggHmumu115_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV115.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-115.log & 
nice ./analyzer ggHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV125.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-125.log & 
nice ./analyzer ggHmumu135_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV135.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-135.log &
nice ./analyzer ggHmumu150_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV150.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-150.log 

nice ./analyzer vbfHmumu115_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV115.root -r 7TeV $OPTIONS >& log_VBF7TeV-115.log &
nice ./analyzer vbfHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV125.root -r 7TeV $OPTIONS >& log_VBF7TeV-125.log &
nice ./analyzer vbfHmumu135_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV135.root -r 7TeV $OPTIONS >& log_VBF7TeV-135.log &
nice ./analyzer vbfHmumu150_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV150.root -r 7TeV $OPTIONS >& log_VBF7TeV-150.log 

wait

### Varied samples
echo "Running Variation Samples..."

for i in /data/uftrig01b/jhugon/hmumu/samples/systematics/*.root; do
  tmp=$(basename $i | sed "s/.*\([0-9]TeV\).*/\1/")
  echo "nice ./analyzer $(basename $i) $i -r $tmp >& log_$(basename $i)"
  nice ./analyzer $(basename $i) $i -r $tmp >& log_$(basename $i)
done

echo "Done."
