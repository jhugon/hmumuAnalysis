#!/bin/bash

nice scons -j4

#TRAININGTREES="true"
#TRAIN="true"

DIR=/data/uftrig01b/digiovan/root/higgs/CMSSW_5_3_3_patch3/V00-01-01/

echo "#######################"
echo "   Running Analyzer"
echo "#######################"
echo `date`

if [ "$TRAININGTREES" = "true" ]; then
echo "#######################"
echo "Creating Training Trees"
echo "#######################"

nice ./analyzer DYJetsToLL.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root --trainingTree backgroundTreeDY.root >& log2 &
nice ./analyzer ttbar.root $DIR/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/TTJets_minimal.root --trainingTree backgroundTreeTT.root >& log2 &

nice ./analyzer ggHmumu125.root $DIR/NtuplesMCPrivateSignal/ggHmumu8TeV125.root --trainingTree signalTreeGG.root
nice ./analyzer vbfHmumu125.root $DIR/NtuplesMCPrivateSignal/vbfHmumu8TeV125.root --trainingTree signalTreeVBF.root
nice ./analyzer zHmumu125.root $DIR/NtuplesMCPrivateSignal/zHmumu8TeV125.root --trainingTree signalTreeZH.root
nice ./analyzer wHmumu125.root $DIR/NtuplesMCPrivateSignal/wHmumu8TeV125.root --trainingTree signalTreeWH.root

nice ./analyzer WW.root $DIR/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root --trainingTree backgroundTreeWW.root
nice ./analyzer WZ.root $DIR/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root --trainingTree backgroundTreeWZ.root
nice ./analyzer ZZ.root $DIR/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root --trainingTree backgroundTreeZZ.root

#nice ./analyzer DYToTauTau.root $DIR/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root --trainingTree backgroundTreeDYToTauTau.root

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
nice ./mvaTrain inclusive.cfg >& logMVAInc
echo "training VBF..."
nice ./mvaTrain vbf.cfg >& logMVAVBF
echo "done training."
echo "#######################"
echo "#######################"
fi

# Run with full MVA
nice ./analyzer DYJetsToLL.root $DIR/NtuplesMCDYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/DYToMuMu_minimal.root >& log2 &
nice ./analyzer ttbar.root $DIR/NtuplesMCTTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/minimal/TTJets_minimal.root >& log2 &

nice ./analyzer ggHmumu125.root $DIR/NtuplesMCPrivateSignal/ggHmumu8TeV125.root
nice ./analyzer vbfHmumu125.root $DIR/NtuplesMCPrivateSignal/vbfHmumu8TeV125.root
nice ./analyzer zHmumu125.root $DIR/NtuplesMCPrivateSignal/zHmumu8TeV125.root
nice ./analyzer wHmumu125.root $DIR/NtuplesMCPrivateSignal/wHmumu8TeV125.root

nice ./analyzer DYToTauTau.root $DIR/NtuplesMCDYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root
nice ./analyzer WW.root $DIR/NtuplesMCWW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root
nice ./analyzer WZ.root $DIR/NtuplesMCWZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root 
nice ./analyzer ZZ.root $DIR/NtuplesMCZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*.root 
#nice ./analyzer WJetsToLNu.root $DIR/NtuplesMCWJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/*.root >& log2 &
#nice ./analyzer QCD.root $DIR/NtuplesMCQCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/*.root >& log2 &

nice ./analyzer SingleMuRun2012Av1.root $DIR/NtuplesDataSingleMuRun2012A-13Jul2012-v1/minimal/SingleMuRun2012A-13Jul2012-v1_minimal.root
#nice ./analyzer SingleMuRun2012Bv1.root $DIR/NtuplesDataSingleMuRun2012B-13Jul2012-v1/*.root
#nice ./analyzer SingleMuRun2012Cv1.root $DIR/NtuplesDataSingleMuRun2012C-PromptReco-v1/*.root
#nice ./analyzer SingleMuRun2012Cv2.root $DIR/NtuplesDataSingleMuRun2012C-PromptReco-v2/*.root

#nice ./analyzer DoubleMuRun2012Av1.root $DIR/NtuplesDataDoubleMuRun2012A-13Jul2012-v1/*.root
#nice ./analyzer DoubleMuRun2012Bv1.root $DIR/NtuplesDataDoubleMuRun2012B-13Jul2012-v4/*.root
#nice ./analyzer DoubleMuRun2012Cv1.root $DIR/NtuplesDataDoubleMuRun2012C-PromptReco-v1/*.root
#nice ./analyzer DoubleMuRun2012Cv2.root $DIR/NtuplesDataDoubleMuRun2012C-PromptReco-v2/*.root

wait

echo "#######################"
echo "#######################"
echo `date`
echo "Done."
