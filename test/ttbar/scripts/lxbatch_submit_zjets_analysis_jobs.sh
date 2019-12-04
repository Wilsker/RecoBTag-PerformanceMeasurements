bsub -q 1nd -J zjetsanalyser_diboson < lxbatch_run_zjets_analysis_Diboson_jobs.py
bsub -q 1nd -J zjetsanalyser_ST < lxbatch_run_zjets_analysis_ST_jobs.py
bsub -q 1nd -J zjetsanalyser_DYHM < lxbatch_run_zjets_analysis_DYHighMass_jobs.py
bsub -q 1nd -J zjetsanalyser_tt < lxbatch_run_zjets_analysis_ttbar_jobs.py
bsub -q 1nd -J zjetsanalyser_DYLM < lxbatch_run_zjets_analysis_DYLowMass_jobs.py
bsub -q 1nd -J zjetsanalyser_dblEG < lxbatch_run_zjets_analysis_Data_DoubleEG_jobs.py
bsub -q 1nd -J zjetsanalyser_dblMu < lxbatch_run_zjets_analysis_Data_DoubleMuon_jobs.py
