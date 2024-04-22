########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file can be used to run the core function with and without log-likelihood 
# calculation and comparison with previously stored results as unit-test.
#
# Copyright 2024, SIMID                                       
########################################################################### #

rm(list=ls())

# load functions and data ----
source('R/main_vaccination.R')
source('R/lib_unit_testing.R')

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# set parameter file
chains_param_file <- 'data/config/MCMCmulti_20230120_d730_belgium_ve00_customQ.csv'

# boolean to enable the comparison of the new output with the previously stored output (FALSE allows to store reference output)
bool_compare   <- TRUE

# run benchmark?
bool_benchmark <- FALSE

# load file
chains_param      = read.table(chains_param_file, sep = ",", header = T)

# select one parameter set
parms      <- unlist(chains_param[1,])
length(parms)

# store parameter names
parms_names = names(parms)

# set rng seed
set.seed(666)

# set plot configuration
plots_code <- TRUE

# set model output configuration
bool_full_output <- TRUE

# set vaccine uptake scheme
vacc_schedule      <- read.table("data/uptake/vaccine_uptake_booster_vSCENmay26_belgium.csv", sep=',',header = T)

# reference data
be_ref_data      <- read.csv('data/covid19_reference_data_20230526_1030_belgium.csv')

# contact data
contact_data <- readRDS(file='data/social_contact_data.rds')


################################################################ #
## RUN ----
################################################################ #

# Deterministic (mean) model ----
print(paste(Sys.time(),"START DETERMINISTIC MODEL (v.x)")); ptm = proc.time()
fit_mean_scm = log_likelihood_model(parms, 
                                contact_data = contact_data,
                                method = "mean",
                                plots = plots_code,
                                parms_names = names(parms),
                                vaccine_uptake = vacc_schedule,
                                be_ref_data = be_ref_data)
print(proc.time()-ptm) ;
if(bool_compare) compare_output(fit_mean_scm,'mean',prev_file_names)

# Stochastic model ----
print(paste(Sys.time(),"START STOCHASTIC MODEL (v.x)")); ptm = proc.time()
fit_stochastic_scm = log_likelihood_model(parms, 
                                      contact_data = contact_data,
                                      method = "stochastic",
                                      plots = plots_code,
                                      parms_names = names(parms),
                                      vaccine_uptake = vacc_schedule,
                                      be_ref_data = be_ref_data)
print(proc.time()-ptm)
if(bool_compare) compare_output(fit_stochastic_scm,'stochastic',prev_file_names)

## MODEL EXPLORATION ----
print(paste(Sys.time(),"EXPLORE STOCHASTIC MODEL (v.x)")); ptm = proc.time()
scm_out = run_model(parms,
                contact_data = contact_data,
                method = "stochastic",
                parms_names = names(parms),
                ndays_sim = 850,
                vaccine_uptake = vacc_schedule,
                be_ref_data = be_ref_data)

print(proc.time()-ptm) ;
if(bool_compare) compare_output(scm_out,'scm_out',prev_file_names)

print(paste(Sys.time(),"EXPLORE DETERMINISTIC MODEL (v.x)")); ptm = proc.time()
scm_mean = run_model(parms,
                contact_data = contact_data,
                method = "mean",
                parms_names = names(parms),
                ndays_sim = 900,
                vaccine_uptake = vacc_schedule,
                be_ref_data = be_ref_data)

print(proc.time()-ptm) ;
if(bool_compare) compare_output(scm_mean,'scm_mean',prev_file_names)

## PRIOR ----
print(paste(Sys.time(),"EXPLORE LOG PRIOR (v.x)")); ptm = proc.time()
scm_prior = log_prior_model(parms,
                            parms_names = names(parms))

if(any(grepl('invalid_param',names(scm_prior)))){
  print(scm_prior$invalid_param)
  scm_prior <- scm_prior[1:2]
}

print(proc.time()-ptm) ;
if(bool_compare) compare_output(scm_prior,'scm_prior',prev_file_names)

## EXTENDED TIME HORIZON
print(paste(Sys.time(),"EXPLORE EXTENDED TIME HORIZON, DETERMINISTIC (v.x)")); ptm = proc.time()
fit_extended_scm = log_likelihood_model(parms, 
                                  contact_data = contact_data,
                                  method = "mean",
                                  plots = TRUE,
                                  parms_names = names(parms),
                                  ndays_sim = 850,
                                  vaccine_uptake = vacc_schedule,
                                  be_ref_data = be_ref_data)
print(proc.time()-ptm) ;
if(bool_compare) compare_output(fit_extended_scm,'extended',prev_file_names)

## ADJUSTED BEHAVOIR
print(paste(Sys.time(),"EXPLORE ADJUSTED BEHAVIOR, DETERMINISTIC (v.x)")); ptm = proc.time()
db_cnt_adjust <- data.frame(value=1.7,
                            date = "2021-03-15",
                            reference_date = "2021-03-02",
                            is_parmanent = TRUE)

scm_behaviour = run_model(parms,
                     contact_data = contact_data,
                     method = "mean",
                     parms_names = names(parms),
                     ndays_sim = 400,
                     vaccine_uptake = vacc_schedule,
                     be_ref_data = be_ref_data,
                     db_cnt_adjust = db_cnt_adjust)

print(proc.time()-ptm) ;
if(bool_compare) compare_output(scm_behaviour,'scm_behaviour',prev_file_names)


################################################################ #
## BENCHMARK ----
################################################################ #

if(bool_benchmark){
# note: reformatting the uptake once, instead of every iteration, saves 30% time! 
# decreased from 1.2s to 0.8s per run on Nov 30, 2021
# increased to 0.94s with waning immunity (default reference data)
# 1.3s with waning immunity and latest reference data (December 1st)
ndays_sim     <- nrow(be_ref_data)
V_mat         <- get_vaccine_protection_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_A_', parms = parms, ndays_sim)
V_mat_dose2   <- get_vaccine_protection_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_B_', parms = parms, ndays_sim)
V_mat_booster <- get_vaccine_protection_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_E_', parms = parms, ndays_sim)

print(paste(Sys.time(),"BENCHMARK STOCHASTIC MODEL (v.x)")); ptm = proc.time()
nrun <- 10
for(i in 1:nrun){
  fit_mean_bench = log_likelihood_model(parms, 
                                     contact_data = contact_data,
                                     method = "mean",
                                     plots = FALSE,
                                     parms_names = names(parms),
                                     be_ref_data = be_ref_data,
                                     V_mat = V_mat,
                                     V_mat_dose2 = V_mat_dose2,
                                     V_mat_booster = V_mat_booster)
}
print((total_time = proc.time()-ptm));
print(total_time/nrun) ;
} else{
  print(paste(Sys.time(),"NO BENCHMARKING THIS TIME"))
}

print(paste(Sys.time(),"WORKBENCH FINISHED"))


################################################################ #
# # reset reference values (optional)
################################################################ #
#rrv()


