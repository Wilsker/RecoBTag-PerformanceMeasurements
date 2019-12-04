#!/bin/bash

INPUTFILE="UltraLegacy_2018Data/combined_files/MC13TeV_TTJets_training_combined.root"

a=(0 1 2)
for i in ${a[@]}; do
    root -b -q "${CMSSW_BASE}/src/RecoBTag/PerformanceMeasurements/test/ttbar/KIN_trainClassifier.C+(\"BDT\",\"${INPUTFILE}\",${i})";
    mv TMVA.root TMVA_${i}.root;
done
