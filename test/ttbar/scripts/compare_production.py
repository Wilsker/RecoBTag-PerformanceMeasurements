#############################
#   compare_production.py
#############################
#  Joshuha Thomas-Wilsker
#   CERN, IHEP Beijing
#############################
# Short python script that
# compares production files
# sizes stored on eos and
# gives warnings if it finds
# any differences.
#############################


import os

temp_eos_ls_cmd_1 = 'eos ls -l /store/group/phys_btag/Commissioning/TTbar/KinMethod_StructuredDir_10Feb2017/crab_Data13TeV_MuonEG_2017D_Nov17/'
temp_os_output_1 = os.popen(temp_eos_ls_cmd_1).read()
file_list_1 = temp_os_output_1.splitlines()

temp_eos_ls_cmd_2 = 'eos ls -l /store/group/phys_btag/Commissioning/TTbar/Jan312017withPU/7c84d07/MuonEG/crab_Data13TeV_MuonEG_2017D_Nov17/180131_021316/0000/'
temp_os_output_2 = os.popen(temp_eos_ls_cmd_2).read()
file_list_2 = temp_os_output_2.splitlines()

logfile1_filesizes = []
logfile2_filesizes = []
if len(file_list_1) != len(file_list_2):
    print 'WARNING! # files in directories is different. Something is amiss!'

for x in xrange(len(file_list_1)):
    file_info_1 = file_list_1[x]
    file_info_2 = file_list_2[x]
    start_of_target_string_1 = file_info_1[39:]
    full_target_string_1 = start_of_target_string_1[:7]
    start_of_target_string_2 = file_info_2[39:]
    full_target_string_2 = start_of_target_string_2[:7]
    if full_target_string_1 != full_target_string_2:
        print 'WARNING! files of different sizes:'
        print 'File 1: ', file_info_1
        print 'file 1 size: ', full_target_string_1
        print 'File 2: ', file_info_2
        print 'file 2 size: ', full_target_string_2
    logfile1_filesizes.append(full_target_string_1)
    logfile2_filesizes.append(full_target_string_2)

#print logfile1_filesizes
#print logfile2_filesizes
