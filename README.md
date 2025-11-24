This repository contains code and data for the paper: Estimating pre-symptomatic transmission potential of influenza A and B viruses in household transmission studies

## Data ##
infA_hk_data.csv: Household data of influenza A\n
infB_hk_data.csv: Household data of influenza B\n
infA_incidence: Influenza A activity in community during the study period

## Model C++ code ##
infA_model.cpp: Influenza A main model\n
infB_model.cpp: Influenza B main model\n
covidhk_model.cpp: SARS-CoV-2 main model

## R scripts to perform model inference ##
infA_main.R\n
infB_main.R\n
covidhk_main.R

## Helper functions ##
helper_funcs.R

## Results ##
30,000 MCMC results: tt_infA.RData tt_infB.RData tt_covidhk.RData\n
Estimated values of model parameters: mcmc_summary_infA.csv mcmc_summary_infB.csv mcmc_summary_covidhk.csv\n
Model adequacy check: adeq_infA.csv adeq_infB.csv adeq_covidhk.csv
