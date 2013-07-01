#!/bin/bash

nice scons -j4

#TRAININGTREES="true"
#TRAIN="true"
# Makes each analyzer only run on 1k events:
#OPTIONS=" -m 1000"
# Makes each analyzer only run on the data/MC comparison range in m_mumu:
#OPTIONS=$OPTIONS" -d"

#SIGNALEFF="true"

DIR=/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/

echo "#######################"
echo "   Running Analyzer"
echo "#######################"
echo `date`

if [ "$TRAININGTREES" = "true" ]; then
echo "#######################"
echo "Creating Training Trees"
echo "#######################"

nice ./analyzer DYJetsToLL_8TeV.root $DIR/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYJetsToLL_minimal.root --trainingTree backgroundTreeDY_8TeV.root -r 8TeV $OPTIONS >& log_DYJets < /dev/null &
#nice ./analyzer DYToMuMu_8TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root --trainingTree backgroundTreeDY_8TeV.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer ttbar_8TeV.root $DIR/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1/minimal/TTJets_minimal.root --trainingTree backgroundTreeTT_8TeV.root -r 8TeV $OPTIONS  >& log_TTJets < /dev/null&

nice ./analyzer ggHmumu125_8TeV.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCGluGlu_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-125_forTraning.root --trainingTree signalTreeGG_8TeV.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer vbfHmumu125_8TeV.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCVBF_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-125_forTraining.root --trainingTree signalTreeVBF_8TeV.root -r 8TeV $OPTIONS 

#nice ./analyzer WW_8TeV.root $DIR/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WW_minimal.root --trainingTree backgroundTreeWW_8TeV.root -r 8TeV $OPTIONS  >& log2 &
#nice ./analyzer WZ_8TeV.root $DIR/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WZ_minimal.root --trainingTree backgroundTreeWZ_8TeV.root -r 8TeV $OPTIONS 
#nice ./analyzer ZZ_9TeV.root $DIR/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/ZZ_minimal.root --trainingTree backgroundTreeZZ_8TeV.root -r 8TeV $OPTIONS 

##nice ./analyzer DYToTauTau_8TeV.root $DIR/NtuplesMCDYToTauTau_M_20_TuneZ2star_8TeV_pythia6_v2_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToTauTau_minimal.root --trainingTree backgroundTreeDYToTauTau_8TeV.root -r 8TeV $OPTIONS

wait

echo "#######################"
echo "#######################"
fi

if [ "$TRAIN" = "true" ]; then
echo "#######################"
echo "    Training MVAs"
echo "#######################"
echo "#######################" >& log2
echo "    Training MVAs" >& log2
echo "#######################" >& log2
echo "training Inclusive..."
nice ./mvaTrain inclusive_8TeV.cfg >& logMVAInc_8TeV
echo "training VBF..."
nice ./mvaTrain vbf_8TeV.cfg >& logMVAVBF_8TeV
echo "done training."
echo "#######################"
echo "#######################"
fi

# Run with full MVA
nice ./analyzer DYJetsToLL_8TeV.root $DIR/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYJetsToLL_minimal.root -r 8TeV $OPTIONS 2>&1 >> log2 &
#nice ./analyzer DYToMuMu_8TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer ttbar_8TeV.root $DIR/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7C-v1/minimal/TTJets_minimal.root -r 8TeV $OPTIONS  2>&1 >> log2 &

nice ./analyzer ggHmumu125_8TeV.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCGluGlu_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-125_forXcheck.root -r 8TeV $OPTIONS  2>&1 >> log2 &
nice ./analyzer vbfHmumu125_8TeV.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCVBF_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-125_forXcheck.root -r 8TeV $OPTIONS 
nice ./analyzer wHmumu125_8TeV.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCPrivateSignal/wHmumu8TeV125.root -r 8TeV $OPTIONS  2>&1 >> log2 &
nice ./analyzer zHmumu125_8TeV.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_5/V00-01-10/NtuplesMCPrivateSignal/zHmumu8TeV125.root -r 8TeV $OPTIONS 

#nice ./analyzer DYToTauTau_8TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root -r 8TeV $OPTIONS  >& log2 &
nice ./analyzer WW_8TeV.root $DIR/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WW_minimal.root -r 8TeV $OPTIONS 
nice ./analyzer WZ_8TeV.root $DIR/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/WZ_minimal.root -r 8TeV $OPTIONS   2>&1 >> log2 &
nice ./analyzer ZZ_8TeV.root $DIR/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/ZZ_minimal.root -r 8TeV $OPTIONS  
#nice ./analyzer WJetsToLNu_8TeV.root $DIR/NtuplesMCWJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/*.root -r 8TeV $OPTIONS >& log2 &
#nice ./analyzer QCD_8TeV.root $DIR/NtuplesMCQCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/minimal/QCD_Pt_20_MuEnrichedPt_15_minimal.root -r 8TeV $OPTIONS >& log2 &

nice ./analyzer SingleMuRun2012Av1.root $DIR/NtuplesDataSingleMuRun2012A-13Jul2012-v1/minimal/SingleMuRun2012A-13Jul2012-v1_minimal.root -r 8TeV $OPTIONS  2>&1 >> log2 &
nice ./analyzer SingleMuRun2012Av1Recover.root $DIR/NtuplesDataSingleMuRun2012A-recover-06Aug2012-v1/minimal/SingleMuRun2012A-recover-06Aug2012-v1_minimal.root -r 8TeV $OPTIONS 2>&1 >> log2 &
nice ./analyzer SingleMuRun2012Bv1.root $DIR/NtuplesDataSingleMuRun2012B-13Jul2012-v1/minimal/SingleMuRun2012B-13Jul2012-v1_minimal.root -r 8TeV $OPTIONS 2>&1 >> log2 &
nice ./analyzer SingleMuRun2012Cv1.root $DIR/NtuplesDataSingleMuRun2012C-24Aug2012-v1/minimal/SingleMuRun2012C-24Aug2012-v1_minimal.root -r 8TeV $OPTIONS
nice ./analyzer SingleMuRun2012Cv2.root $DIR/NtuplesDataSingleMuRun2012C-PromptReco-v2/minimal/SingleMuRun2012C-PromptReco-v2_minimal.root -r 8TeV $OPTIONS  2>&1 >> log2 &
nice ./analyzer SingleMuRun2012D.root $DIR/NtuplesDataSingleMuRun2012D-PromptReco-v1/minimal/SingleMuRun2012D-PromptReco-v1_minimal.root -r 8TeV $OPTIONS 

#nice ./analyzer SingleMuRun2012Dv1-22Jan2013.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_9/V00-01-13/NtuplesDataSingleMuRun2012D-22Jan2013-v1/minimal/SingleMuRun2012D-13Jan2013-v1_minimal.root -r 8TeV $OPTIONS  >& log2 &
#nice ./analyzer SingleMuRun2012Av1-22Jan2013.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_9/V00-01-13/NtuplesDataSingleMuRun2012A-22Jan2013-v1/minimal/SingleMuRun2012A-13Jan2013-v1_minimal.root -r 8TeV $OPTIONS 
#nice ./analyzer SingleMuRun2012Bv1-22Jan2013.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_9/V00-01-13/NtuplesDataSingleMuRun2012B-22Jan2013-v1/minimal/SingleMuRun2012B-13Jan2013-v1_minimal.root -r 8TeV $OPTIONS 
#nice ./analyzer SingleMuRun2012Cv1-22Jan2013.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_9/V00-01-13/NtuplesDataSingleMuRun2012C-22Jan2013-v1/minimal/SingleMuRun2012C-13Jan2013-v1_minimal.root -r 8TeV $OPTIONS 
#
#wait

wait

######################################
### for signal efficiency calculation
######################################

if [ "$SIGNALEFF" = "true" ]; then

## Official GluGlu Higgs to MuMu
nice ./analyzer ggHmumu115_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-115_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-115.root -r 8TeV $OPTIONS >& log_GluGlu-115.log & 
nice ./analyzer ggHmumu120_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-120_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-120.root -r 8TeV $OPTIONS >& log_GluGlu-120.log & 
nice ./analyzer ggHmumu124_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-124_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-124.root -r 8TeV $OPTIONS >& log_GluGlu-124.log & 
nice ./analyzer ggHmumu124p5_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-124p5_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-124p5.root -r 8TeV $OPTIONS >& log_GluGlu-124p5.log  
nice ./analyzer ggHmumu125_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-125_forXcheck.root -r 8TeV $OPTIONS >& log_GluGlu-125.log 

nice ./analyzer ggHmumu125p5_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-125p5_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-125p5.root -r 8TeV $OPTIONS >& log_GluGlu-125p5.log & 
nice ./analyzer ggHmumu126_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-126_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-126.root -r 8TeV $OPTIONS >& log_GluGlu-126.log & 
nice ./analyzer ggHmumu126p5_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-126p5_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-126p5.root -r 8TeV $OPTIONS >& log_GluGlu-126p5.log  

nice ./analyzer ggHmumu127_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-127_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-127.root -r 8TeV $OPTIONS >& log_GluGlu-127.log & 
nice ./analyzer ggHmumu128_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-128_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-128.root -r 8TeV $OPTIONS >& log_GluGlu-128.log & 
nice ./analyzer ggHmumu130_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-130_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-130.root -r 8TeV $OPTIONS >& log_GluGlu-130.log & 
nice ./analyzer ggHmumu135_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-135_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-135.root -r 8TeV $OPTIONS >& log_GluGlu-135.log  
###nice ./analyzer ggHmumu140_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-140_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-140.root -r 8TeV $OPTIONS >& log_GluGlu-140.log & 
nice ./analyzer ggHmumu145_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-145_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-145.root -r 8TeV $OPTIONS >& log_GluGlu-145.log & 
nice ./analyzer ggHmumu150_8TeV.root $DIR/NtuplesMCGluGlu_HToMM_M-150_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/GluGlu_HToMM_M-150.root -r 8TeV $OPTIONS >& log_GluGlu-150.log & 
nice ./analyzer ggHmumu155_8TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu8TeV155.root -r 8TeV $OPTIONS >& log_GluGlu-155.log & 
nice ./analyzer ggHmumu160_8TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu8TeV160.root -r 8TeV $OPTIONS >& log_GluGlu-160.log  

## Official VBF Higgs to MuMu
nice ./analyzer vbfHmumu115_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-115_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-115.root -r 8TeV $OPTIONS >& log_VBF-115.log & 
nice ./analyzer vbfHmumu120_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-120_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-120.root -r 8TeV $OPTIONS >& log_VBF-120.log & 
nice ./analyzer vbfHmumu124_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-124_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-124.root -r 8TeV $OPTIONS >& log_VBF-124.log & 
####nice ./analyzer vbfHmumu124p5_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-124p5_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-124p5.root -r 8TeV $OPTIONS >& log_VBF-124p5.log & 
nice ./analyzer vbfHmumu125_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-125_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-125_forXcheck.root -r 8TeV $OPTIONS >& log_VBF-125.log
nice ./analyzer vbfHmumu125p5_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-125p5_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-125p5.root -r 8TeV $OPTIONS >& log_VBF-125p5.log & 
nice ./analyzer vbfHmumu126_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-126_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-126.root -r 8TeV $OPTIONS >& log_VBF-126.log & 
nice ./analyzer vbfHmumu126p5_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-126p5_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-126p5.root -r 8TeV $OPTIONS >& log_VBF-126p5.log 

nice ./analyzer vbfHmumu127_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-127_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-127.root -r 8TeV $OPTIONS >& log_VBF-127.log & 
nice ./analyzer vbfHmumu128_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-128_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-128.root -r 8TeV $OPTIONS >& log_VBF-128.log & 
###nice ./analyzer vbfHmumu130_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-130_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-130.root -r 8TeV $OPTIONS >& log_VBF-130.log & 
nice ./analyzer vbfHmumu135_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-135_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-135.root -r 8TeV $OPTIONS >& log_VBF-135.log 
nice ./analyzer vbfHmumu140_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-140_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-140.root -r 8TeV $OPTIONS >& log_VBF-140.log & 
nice ./analyzer vbfHmumu145_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-145_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-145.root -r 8TeV $OPTIONS >& log_VBF-145.log & 
nice ./analyzer vbfHmumu150_8TeV.root $DIR/NtuplesMCVBF_HToMM_M-150_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/VBF_HToMM_M-150.root -r 8TeV $OPTIONS >& log_VBF-150.log 
nice ./analyzer vbfHmumu155_8TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu8TeV155.root -r 8TeV $OPTIONS >& log_VBF-155.log & 
nice ./analyzer vbfHmumu160_8TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu8TeV160.root -r 8TeV $OPTIONS >& log_VBF-160.log & 

wait

fi

tar czf result.tgz ggHmumu*.root vbfHmumu*.root zHmumu*.root wHmumu*.root ttbar*.root DY*.root WW*.root WZ*.root ZZ*.root SingleMu*.root

echo "#######################"
echo "#######################"
echo `date`
echo "Done."
