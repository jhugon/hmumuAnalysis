#!/bin/bash

nice scons -j4

#TRAININGTREES="true"
#TRAIN="true"
# Makes each analyzer only run on 1k events:
#OPTIONS=" -m 1000"
# Makes each analyzer only run on the data/MC comparison range in m_mumu:
#OPTIONS=$OPTIONS" -d"

#SIGNALEFF="true"

DIR=/data/uftrig01b/digiovan/root/higgs/CMSSW_4_4_5/V00-01-10/

echo "#######################"
echo "   Running Analyzer"
echo "#######################"
echo `date`

if [ "$TRAININGTREES" = "true" ]; then
echo "#######################"
echo "Creating Training Trees"
echo "#######################"

nice ./analyzer DYJetsToLL_7TeV.root $DIR/NtuplesMCDYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/DYJetsToLL_minimal.root --trainingTree backgroundTreeDY_7TeV.root -r 7TeV $OPTIONS 2>&1 >> log7TeV2 &
#nice ./analyzer DYToMuMu_7TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START44_V9B-v1/minimal/DYToMuMu_minimal.root --trainingTree backgroundTreeDY_7TeV.root -r 7TeV $OPTIONS  2>&1 >> log7TeV2 &
nice ./analyzer ttbar_7TeV.root $DIR/NtuplesMCTTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/TTJets_minimal.root --trainingTree backgroundTreeTT_7TeV.root -r 7TeV $OPTIONS  2>&1 >> log7TeV2 &

nice ./analyzer ggHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/testForIvan/ggHmumu7TeV125/ggHmmu7TeV125_fortraining_big.root --trainingTree signalTreeGG_7TeV.root -r 7TeV $OPTIONS 
nice ./analyzer vbfHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/testForIvan/vbfHmumu7TeV125/vbfHmmu7TeV125_fortraining_big.root --trainingTree signalTreeVBF_7TeV.root -r 7TeV $OPTIONS 

#nice ./analyzer WW_7TeV.root $DIR/NtuplesMCWW_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WW_minimal.root --trainingTree backgroundTreeWW_7TeV.root -r 7TeV $OPTIONS 
#nice ./analyzer WZ_7TeV.root $DIR/NtuplesMCWZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WZ_minimal.root --trainingTree backgroundTreeWZ_7TeV.root -r 7TeV $OPTIONS 
#nice ./analyzer ZZ_7TeV.root $DIR/NtuplesMCZZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/ZZ_minimal.root --trainingTree backgroundTreeZZ_7TeV.root -r 7TeV $OPTIONS 

#nice ./analyzer DYToTauTau_7TeV.root $DIR/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/DYToTauTau_minimal.root --trainingTree backgroundTreeDYToTauTau_7TeV.root -r 7TeV $OPTIONS

wait

echo "#######################"
echo "#######################"
fi

if [ "$TRAIN" = "true" ]; then
echo "#######################"
echo "    Training MVAs"
echo "#######################"
echo "#######################" 2>&1 >> log7TeV2
echo "    Training MVAs" 2>&1 >> log7TeV2
echo "#######################" 2>&1 >> log7TeV2
echo "training Inclusive..."
nice ./mvaTrain inclusive_7TeV.cfg 2>&1 >> logMVAInc7TeV
echo "training VBF..."
nice ./mvaTrain vbf_7TeV.cfg 2>&1 >> logMVAVBF7TeV
echo "done training."
echo "#######################"
echo "#######################"
fi

# Run with full MVA
nice ./analyzer DYJetsToLL_7TeV.root $DIR/NtuplesMCDYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/DYJetsToLL_minimal.root -r 7TeV $OPTIONS  2>&1 >> log7TeV2 &
#nice ./analyzer DYToMuMu_7TeV.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START44_V9B-v1/minimal/DYToMuMu_minimal.root -r 7TeV $OPTIONS  2>&1 >> log7TeV2 &
nice ./analyzer ttbar_7TeV.root $DIR/NtuplesMCTTJets_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/TTJets_minimal.root -r 7TeV $OPTIONS  2>&1 >> log7TeV2 &

nice ./analyzer ggHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/testForIvan/ggHmumu7TeV125/ggHmmu7TeV125_forxcheck_big.root -r 7TeV $OPTIONS 
nice ./analyzer vbfHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/testForIvan/vbfHmumu7TeV125/vbfHmmu7TeV125_forxcheck_big.root -r 7TeV $OPTIONS 
nice ./analyzer zHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/zHmumu7TeV125.root -r 7TeV $OPTIONS 
nice ./analyzer wHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/wHmumu7TeV125.root -r 7TeV $OPTIONS 

#nice ./analyzer DYToTauTau_7TeV.root $DIR/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/DYToTauTau_minimal.root -r 7TeV $OPTIONS 
nice ./analyzer WW_7TeV.root $DIR/NtuplesMCWW_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WW_minimal.root -r 7TeV $OPTIONS 
nice ./analyzer WZ_7TeV.root $DIR/NtuplesMCWZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WZ_minimal.root -r 7TeV $OPTIONS  
nice ./analyzer ZZ_7TeV.root $DIR/NtuplesMCZZ_TuneZ2_7TeV_pythia6_tauola_Fall11-PU_S6_START44_V9B-v1/minimal/ZZ_minimal.root -r 7TeV $OPTIONS  
#nice ./analyzer WJetsToLNu_7TeV.root $DIR/NtuplesMCWJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall11-PU_S6_START44_V9B-v1/minimal/WJetsToLNu_minimal.root -r 7TeV $OPTIONS 2>&1 >> log7TeV2 &
#nice ./analyzer QCD_7TeV.root $DIR/NtuplesMCQCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6_Fall11-PU_S6_START44_V9B-v1/minimal/QCD_Pt_20_MuEnrichedPt_15_minimal.root -r 7TeV $OPTIONS 2>&1 >> log7TeV2 &

nice ./analyzer ggHmumu125TrainingSample_7TeV.root $DIR/NtuplesMCPrivateSignal/testForIvan/ggHmumu7TeV125/ggHmmu7TeV125_fortraining_big.root -r 7TeV $OPTIONS  2>&1 >> logTrainingSampleVBF27 &
nice ./analyzer vbfHmumu125TrainingSample_7TeV.root $DIR/NtuplesMCPrivateSignal/testForIvan/vbfHmumu7TeV125/vbfHmmu7TeV125_fortraining_big.root -r 7TeV $OPTIONS 2>&1 >> logTrainingSampleVBF7 &

nice ./analyzer SingleMuRun2011Av1.root $DIR/NtuplesDataSingleMuRun2011A-08Nov2011-v1/minimal/SingleMuRun2011A-08Nov2011-v1_minimal.root -r 7TeV $OPTIONS 2>&1 >> log7TeV2 &
nice ./analyzer SingleMuRun2011Bv1.root $DIR/NtuplesDataSingleMuRun2011B-19Nov2011-v1/minimal/SingleMuRun2011B-19Nov2011-v1_minimal.root -r 7TeV $OPTIONS 

wait

######################################
### for signal efficiency calculation
######################################

if [ "$SIGNALEFF" = "true" ]; then

## Private GluGlu Higgs to MuMu
nice ./analyzer ggHmumu115_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV115.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-115.log & 
nice ./analyzer ggHmumu120_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV120.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-120.log & 
nice ./analyzer ggHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/testForIvan/ggHmumu7TeV125/ggHmmu7TeV125_forxcheck_big.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-125.log & 
nice ./analyzer ggHmumu130_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV130.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-130.log &
nice ./analyzer ggHmumu135_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV135.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-135.log 
nice ./analyzer ggHmumu140_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV140.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-140.log & 
nice ./analyzer ggHmumu145_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV145.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-145.log & 
nice ./analyzer ggHmumu150_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV150.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-150.log &
nice ./analyzer ggHmumu155_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV155.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-155.log & 
nice ./analyzer ggHmumu160_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/ggHmmu7TeV160.root -r 7TeV $OPTIONS >& log_GluGlu7TeV-160.log 

nice ./analyzer vbfHmumu115_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV115.root -r 7TeV $OPTIONS >& log_VBF7TeV-115.log & 
nice ./analyzer vbfHmumu120_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV120.root -r 7TeV $OPTIONS >& log_VBF7TeV-120.log & 
nice ./analyzer vbfHmumu125_7TeV.root $DIR/NtuplesMCPrivateSignal/testForIvan/vbfHmumu7TeV125/vbfHmmu7TeV125_forxcheck_big.root -r 7TeV $OPTIONS >& log_VBF7TeV-125.log & 
nice ./analyzer vbfHmumu130_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV130.root -r 7TeV $OPTIONS >& log_VBF7TeV-130.log & 
nice ./analyzer vbfHmumu135_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV135.root -r 7TeV $OPTIONS >& log_VBF7TeV-135.log 
nice ./analyzer vbfHmumu140_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV140.root -r 7TeV $OPTIONS >& log_VBF7TeV-140.log & 
nice ./analyzer vbfHmumu145_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV145.root -r 7TeV $OPTIONS >& log_VBF7TeV-145.log & 
nice ./analyzer vbfHmumu150_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV150.root -r 7TeV $OPTIONS >& log_VBF7TeV-150.log & 
nice ./analyzer vbfHmumu155_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV155.root -r 7TeV $OPTIONS >& log_VBF7TeV-155.log & 
nice ./analyzer vbfHmumu160_7TeV.root $DIR/NtuplesMCPrivateSignal/HPC/100K/vbfHmmu7TeV160.root -r 7TeV $OPTIONS >& log_VBF7TeV-160.log 

wait


fi


tar czf result.tgz ggHmumu*.root vbfHmumu*.root zHmumu*.root wHmumu*.root ttbar*.root DY*.root WW*.root WZ*.root ZZ*.root SingleMu*.root

echo "#######################"
echo "#######################"
echo `date`
echo "Done."
