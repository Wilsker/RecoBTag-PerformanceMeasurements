#!/bin/bash
#Relocate CRAB outputs to improve directory structure.
#=========================================
#CRAB output
lnfpath=/store/group/phys_btag/Commissioning/TTbar/Jan152017withPU/7c84d07/

#New dir
newlfnpath=/store/group/phys_btag/Commissioning/TTbar/KinMethod_StructuredDir_30Jan2017/

#CRAB output dataset titles.
datasets=(
"
TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8
TTToHadronic_TuneCP5_13TeV-powheg-pythia8
TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8
ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8
ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8
WW_TuneCP5_13TeV-pythia8
WZ_TuneCP5_13TeV-pythia8
ZZ_TuneCP5_13TeV-pythia8
DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8
DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8
DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8
DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8
DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8
DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8
DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8
"
#TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8
#TT_TuneCUETP8M2T4_13TeV-powheg-fsrup-pythia8
#TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8
#TT_TuneCUETP8M2T4_13TeV-powheg-isrup-pythia8
#TT_TuneCUETP8M2T4_mtop1695_13TeV-powheg-pythia8
#TT_TuneCUETP8M2T4_mtop1755_13TeV-powheg-pythia8
#TT_TuneCUETP8M2T4down_13TeV-powheg-pythia8
#TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8
)

#Dataset titles for new directory structure.
newdatasets=(
MC13TeV_TTJets_DL
MC13TeV_TTJets_AH
MC13TeV_TTJets_SL
MC13TeV_ST_TW_antitop
MC13TeV_ST_TW_top
MC13TeV_WW
MC13TeV_WZ
MC13TeV_ZZ
MC13TeV_DYJetsToLL-HT_100-200
MC13TeV_DYJetsToLL-HT_200-400
MC13TeV_DYJetsToLL-HT_400-600
MC13TeV_DYJetsToLL-HT_600-800
MC13TeV_DYJetsToLL-HT_800-1200
MC13TeV_DYJetsToLL-HT_2500-inf
MC13TeV_DYJetsToLL-FXFX
#MC13TeV_TTJets_fsrdown
#MC13TeV_TTJets_fsrup
#MC13TeV_TTJets_isrdown
#MC13TeV_TTJets_isrup
#MC13TeV_TTJets_m169v5
#MC13TeV_TTJets_m175v5
#MC13TeV_TTJets_CUETP8M2T4down
#MC13TeV_TTJets_CUETP8M2T4up
)

pos=0
for d in $datasets; do
 dataset=$d
 newdataset=${newdatasets[$pos]}
 #Create new folder
 #eos rm -r $newlfnpath$newdataset
 xrd eoscms mkdir $newlfnpath$newdataset
 echo 'Copying to new directory:'
 echo $newlfnpath$newdataset

 #Get full path
 echo 'Directory:'
 eos ls $lnfpath$dataset
 eos ls $lnfpath$dataset > temp.txt

 newfolder1=`cat temp.txt`
 echo 'Folder:'
 eos ls $lnfpath$dataset/$newfolder1
 eos ls $lnfpath$dataset/$newfolder1 > temp.txt

 newfolder2=`cat temp.txt`
 echo 'Subfolder:'
 eos ls $lnfpath$dataset/$newfolder1/$newfolder2
 eos ls $lnfpath$dataset/$newfolder1/$newfolder2 > temp.txt
 subfolds=`cat temp.txt`

 #Get files
 for sf in $subfolds; do
   echo 'File:'
   echo $lnfpath$dataset/$newfolder1/$newfolder2/$sf
   eos ls $lnfpath$dataset/$newfolder1/$newfolder2/$sf > temp.txt
   files=`cat temp.txt | sed -e "s/.*\///" | grep root | awk '{ print $NF }'`
   #Copy files
   for f in $files; do
     echo 'cp file:'
     echo $lnfpath$dataset/$newfolder1/$newfolder2/$sf/$f
     echo 'into:'
     echo $newlfnpath$newdataset/$f
     xrd eoscms cp $lnfpath$dataset/$newfolder1/$newfolder2/$sf/$f $newlfnpath$newdataset/$f
   done
  done
  let pos=pos+1
 done
