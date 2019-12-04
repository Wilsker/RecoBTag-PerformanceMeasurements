#!/bin/bash
#bsub -q 8nh -J BTV_ttbarAna_ttbarAH < lxbatchSubmit_ttbarAna_ttbarAH.sh
#bsub -q 8nh -J BTV_ttbarAna_ttbarSL < lxbatchSubmit_ttbarAna_ttbarSL.sh
bsub -q 1nd -J BTV_ttbarAna_ttbarDL < lxbatchSubmit_ttbarAna_ttbarDL.sh
