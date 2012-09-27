#!/bin/bash

nice scons -j4

#TRAININGTREES="true"
#TRAIN="true"

echo "#######################"
echo "   Running Analyzer"
echo "#######################"
echo `date`

if [ "$TRAININGTREES" = "true" ]; then
echo "#######################"
echo "Creating Training Trees"
echo "#######################"
nice ./analyzer ggHmumu125.root /data/uftrig01b/jhugon/hmumu/privateSignalV2/ggHmumu125.root --trainingTree signalTreeGG.root
nice ./analyzer vbfHmumu125.root /data/uftrig01b/jhugon/hmumu/privateSignalV2/vbfHmumu125.root --trainingTree signalTreeVBF.root
nice ./analyzer zHmumu125.root /data/uftrig01b/jhugon/hmumu/privateSignalV2/zHmumu125.root --trainingTree signalTreeZH.root
nice ./analyzer wHmumu125.root /data/uftrig01b/jhugon/hmumu/privateSignalV2/wHmumu125.root --trainingTree signalTreeWH.root

nice ./analyzer DYJetsToLL.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root --trainingTree backgroundTreeDY.root
nice ./analyzer DYToTauTau.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root --trainingTree backgroundTreeDYToTauTau.root
nice ./analyzer WZ.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root --trainingTree backgroundTreeWZ.root
nice ./analyzer ZZ.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root --trainingTree backgroundTreeZZ.root
nice ./analyzer ttbar.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_2_6/NtuplesMCTTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/*.root --trainingTree backgroundTreeTT.root

echo "#######################"
echo "#######################"
fi

if [ "$TRAIN" = "true" ]; then
echo "#######################"
echo "    Training MVAs"
echo "#######################"
echo "training Inclusive..."
./mvaTrain inclusive.cfg >& logMVAInc
echo "training VBF..."
./mvaTrain vbf.cfg >& logMVAVBF
echo "done training."
echo "#######################"
echo "#######################"
fi

# Run with full MVA
nice ./analyzer ggHmumu125.root /data/uftrig01b/jhugon/hmumu/privateSignalV2/ggHmumu125.root
nice ./analyzer vbfHmumu125.root /data/uftrig01b/jhugon/hmumu/privateSignalV2/vbfHmumu125.root
nice ./analyzer zHmumu125.root /data/uftrig01b/jhugon/hmumu/privateSignalV2/zHmumu125.root
nice ./analyzer wHmumu125.root /data/uftrig01b/jhugon/hmumu/privateSignalV2/wHmumu125.root

nice ./analyzer DYJetsToLL.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root
nice ./analyzer DYToTauTau.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesMCDYToTauTau_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root
nice ./analyzer WZ.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root
nice ./analyzer ZZ.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root
nice ./analyzer ttbar.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_2_6/NtuplesMCTTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1/*.root

nice ./analyzer SingleMuRun2012Av1.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesDataSingleMuRun2012A-13Jul2012-v1/*.root
nice ./analyzer SingleMuRun2012Bv1.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesDataSingleMuRun2012B-13Jul2012-v1/*.root
#nice ./analyzer SingleMuRun2012Cv1.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesDataSingleMuRun2012C-PromptReco-v1/*.root
#nice ./analyzer SingleMuRun2012Cv2.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesDataSingleMuRun2012C-PromptReco-v2/*.root

#nice ./analyzer DoubleMuRun2012Av1.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesDataDoubleMuRun2012A-13Jul2012-v1/*.root
#nice ./analyzer DoubleMuRun2012Bv1.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesDataDoubleMuRun2012B-13Jul2012-v4/*.root
#nice ./analyzer DoubleMuRun2012Cv1.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesDataDoubleMuRun2012C-PromptReco-v1/*.root
#nice ./analyzer DoubleMuRun2012Cv2.root /data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/NtuplesDataDoubleMuRun2012C-PromptReco-v2/*.root

echo "#######################"
echo "#######################"
echo `date`
echo "Done."
