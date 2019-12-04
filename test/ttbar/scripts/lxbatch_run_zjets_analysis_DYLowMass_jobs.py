#!/bin/bash
# Lxplus Batch Job Script
CMSSW_SRC=/afs/cern.ch/work/j/jthomasw/private/IHEP/BTV/CMSSW/ntuple-production/CMSSW_9_4_1/src
WORKING_DIR=/afs/cern.ch/work/j/jthomasw/private/IHEP/BTV/CMSSW/ntuple-production/CMSSW_9_4_1/src/RecoBTag/PerformanceMeasurements/test/zjets
echo $CMSSW_SRC
export CMSSW_BASE=$CMSSW_SRC
cd $CMSSW_SRC
eval `scramv1 runtime -sh` # Establishes runtime env. of local shell on node.
cd $WORKING_DIR
echo $WORKING_DIR
python runZJetsAnalysis.py -i /store/group/phys_btag/Commissioning/ZJets/KINMistag-2018-06-09_StructuredDir/ -o ${WORKING_DIR}/output_ntuples_2018-08-09/ -j ${WORKING_DIR}/data/samples_Run2018_Mistag_DYLowMass.json
