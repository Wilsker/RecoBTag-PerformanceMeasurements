###########################################################
#           create_lxbatch_submission_scripts.py
###########################################################
#                Joshuha Thomas-Wilsker
#                  CERN, IHEP Beijing
###########################################################
# Short script to create lxbatch submission scripts used in
# bsub_lxbatch_ttbarAna.sh. Requires the script
# lxbatchSubmit_ttbarAna.sh which is used as a template. To
# create more submission scripts, add entries to
# filename_suffix list. Entry should match the suffix of a
# file in data directory:
# e.g. ${WORKING_DIR}/data/samples_Run2017_MC_ttbarAH.json
###########################################################
# Usage: python create_lxbatch_submission_scripts.py
###########################################################

import os, sys, subprocess

filename_suffix = ['ttbarAH','ttbarSL','ttbarDL']
template_file_name = 'lxbatchSubmit_ttbarAna.sh'
input_file = open(template_file_name,'r')
for suff in filename_suffix:
    output_file_name = 'lxbatchSubmit_ttbarAna_'+suff+'.sh'
    output_file = open(output_file_name,'w')
    input_file.seek(0)
    for line in input_file:
        if 'runTTbarAnalysis.py' in line:
            line = line.replace('ttbarAH',suff)
        output_file.write(line)
    output_file.close()
input_file.close()
