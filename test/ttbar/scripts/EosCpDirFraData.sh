#!/bin/bash


#=========================================
#Relocate CRAB outputs to improve directory structure.
#=========================================

#CRAB output
lnfpath=/store/group/phys_btag/Commissioning/TTbar/TTbarRt8023MVA_Final/
#Directory for new structure.
newlfnpath=/store/group/phys_btag/Commissioning/TTbar/JTW_KIN_ROMEO/

datasets=(
#MuonEG_2016B
#MuonEG_2016C
#MuonEG_2016D
#MuonEG_2016E
#MuonEG_2016F
#MuonEG_2016G
MuonEG_2016H
)

pos=0
for d in $datasets; do
  echo $pos
  dataset=${datasets[$pos]}

  #Get full path
  eos ls $lnfpath$dataset*
  eos ls $lnfpath$dataset* > temp.txt

  multicrabsubfolds=`cat temp.txt`

  for rf in $multicrabsubfolds; do
    echo 'cp files from:'
    echo $lnfpath$rf
    echo ' into: '
    echo $newlfnpath/$dataset/$rf

    #Copy files
    xrd eoscms cp $lnfpath/$rf $newlfnpath/$dataset/$rf
  done
let pos=$((pos+1))
done
