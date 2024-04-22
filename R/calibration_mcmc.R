########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file can be used to estimate model parameters using MCMC
#
# This file can be executed with command line arguments. For example:
# - terminal: Rscript R/calibration_mcmc.R 1 out1 1000 10 2020wave1 belgium
# - terminal: Rscript R/calibration_mcmc.R 1 out1 1000 10 2020wave1 belgium &
# - within R: source('R/calibration_mcmc.R')
# - within R: system('Rscript R/calibration_mcmc.R 1 out1 1000 10 2020wave1 belgium &')
#
# Copyright 2024, SIMID                       
########################################################################### #

# Note: when "executing" an R-file, this is called line by line, which is error prone 
# when one edits the file while running a simulation (in the background). To enable 
# multiple runs of this script the same time, the code is captured in a function and 
# this function is called at the end of this script with the command line arguments.

# parse command line arguments
args = commandArgs(trailingOnly=TRUE)

# clear workspace (except the command line arguments)
rm(list=ls()[ls()!='args'])

# load functions and data ----
source('R/main_vaccination.R')
source('R/lib_calibration.R')

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# main function parameters:
# - f_task_id           the id of the current task/run to conduct one MCMC chain (0:n)
# - f_subdir            the name of the subdirectory to store the model output
# - f_n_iter            the number of MCMC iterations 
# - f_n_period          the MCMC period
# - f_param             the parameter selection
# - f_region            the region to calibrate the model for (wrt previous config files and reference data)
# - f_filename_pattern  to grasp the configuration file with initial values

# uncomment for debugging:
#f_task_id = 1; f_subdir=NA; f_n_iter = 10; f_n_period = 10; f_param = '2020q1'; f_region = "belgium"; f_filename_pattern = NA;f_social_contact_data = NA; f_filename_pattern = NA
run <- function(f_task_id = 1,
                f_subdir='SCM_V28', 
                f_n_iter = 10, 
                f_n_period = 10, 
                f_param = '2020wave1',
                f_region = "belgium",
                f_ref_data_file_name = NA,
                f_social_contact_data = NA,
                f_filename_pattern = NA){
  
  
  
# set default run tag and output dir           
output_tag <- paste0('MCMC_',f_subdir)

# set start parameter filename pattern
filename_pattern <- 'data/config/MCMCmulti_20230120_d730_belgium_ve00_customQ' # belgium

# set directory with vaccine uptake files (region-specific selection is done later)
vaccine_uptake_dir     <- 'data/uptake/uptake_vac_voc'
vaccine_uptake_pattern <- 'vSCENjan09_belgium'

# set default MCMC parameters
n_iter   <- as.numeric(f_n_iter)
n_period <- as.numeric(f_n_period)
n_repeat <- 10
n_status <- 1
sel_crit <- 'crit2'   # the log-likelihood criteria from the stochastic model
s_deviance <- 0.1

# set boolean to plot model output (using initial and final parameters)
bool_plot_model_output <- TRUE

################################################################ #
# LOAD AND PRE-PROCESS ----
################################################################ #
f_chain_id <- as.numeric(f_task_id) + is_vsc_vaughan_cluster() # adjusted for TORQUE 

run_tag <- paste(f_region, f_param, f_chain_id,sep='_')

# create small delay when running different chains
if(as.numeric(f_chain_id)){
  Sys.sleep(f_chain_id)
}

# optional: adjust filname pattern for the parameter file
if(exists('f_filename_pattern') && !is.na(f_filename_pattern)){                          
  filename_pattern   <- f_filename_pattern
} 

# reference data
be_ref_data <- get_latest_incidence_data(sel_region = f_region)
if(!is.na(f_ref_data_file_name)){
  be_ref_data      <- read.table(f_ref_data_file_name,sep=',',header=T)
  be_ref_data$date <- as.Date(be_ref_data$date)
}

# social contact data
contact_data <- readRDS('data/social_contact_data.rds')
if(!is.na(f_social_contact_data)){
  contact_data  <- readRDS(paste0('data/', f_social_contact_data, '.rds'))
}

# set default time horizon
ndays_sim <- sim_date2day(max(be_ref_data$date[!is.na(be_ref_data$hospital_admissions_other)]))
  
# option to use command line MCMC paramters
if(exists('f_n_iter') && !is.na(f_n_iter)){ n_iter <- as.numeric(f_n_iter) }  
if(exists('f_n_period') && !is.na(f_n_period)){ n_period <- as.numeric(f_n_period) }  

# get parameter file(s) to start from (selection is done later)
chains_param_files = dir('data/config/',full.names = T,recursive = T) 

if(is_vsc_vaughan_cluster()){
  param_file_dirs = dir(get_output_dir(),full.names = T,recursive = F,pattern = paste0("MCMC_",filename_pattern,".*_final")) 
  param_file_dirs <- param_file_dirs[!grepl('final_pdf',param_file_dirs)]
  output_param_files = dir(param_file_dirs,full.names = T,recursive = T,pattern = "MCMCmulti")
  chains_param_files <- c(chains_param_files,output_param_files)
}

# option to select region-specific file
chains_param_file <- ifelse(any(grepl(f_region,chains_param_files)),
                            chains_param_files[grepl(f_region,chains_param_files)],
                            chains_param_files[1])

# other option to select a specific parameter file
chains_param_file <- ifelse(any(grepl(filename_pattern,chains_param_files)),
                            chains_param_files[grepl(filename_pattern,chains_param_files)],
                            chains_param_files[1])

# read parameter values
parms_chains      = read.table(chains_param_file, sep = ",", header = T)
print(chains_param_file)
dim(parms_chains)

# select parameters set (option to vary the parameter set using the chain_id)
parms      <- unlist(parms_chains[min(as.numeric(f_chain_id),nrow(parms_chains)),])

# modify q parameters?
parms <- adjust_stage_param(parms,contact_data)

# set parameter names
parms_names = names(parms)

# option to aggregate contact stages
opt_aggr_q_stages <- NA

# update region_id in param
parms['region_id'] <- get_region_id(f_region)

# check names
if(any(is.na(names(parms)))){
  stop('ISSUE WITH PARAMETER NAMES... PLEASE CHECK FOR "NA" ')
}

################################################################ #
# SELECT MODEL PARAMETERS TO ESTIMATE ----
################################################################ #

# DEFAULT PARAMETER SET: all (except "ndays_calibration" and "h, "region_id" and mcmc meta info)
parms_names_estim <- parms_names[-which(parms_names %in% c("ndays_calibration","h","region_id") | grepl('mcmc_',parms_names))]
num_param_all     <- length(parms_names_estim)
num_stages        <- min(5,nrow(contact_data$db_C_sim))

# option: 2020wave1  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2020wave1'){
  ndays_sim <- 122 
  
  # check and limit adjusted omega parameters
  flag_omega <- grepl('log_omega',parms_names)
  parms[flag_omega & parms < log(1/7)] <- log(1/2)
  parms[flag_omega & parms > log(2)]   <- log(1/2)
  parms_chains[,flag_omega]            <- log(1/2)
  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    sel_id = 1:min(which(grepl('log_n0',parms_names))),
                                    tag_list = paste0('stage',1:5,'_'),
                                    ignore_list = c('log_phi1',
                                                    'log_beta1',
                                                    'log_delta3',
                                                    'log_delta4',
                                                    'log_mu')
  )
  sel_crit <- 'crit1'
}

# option: 2020wave2----
if(exists('f_param') && !is.na(f_param) && f_param == '2020wave2'){
  ndays_sim <- 306 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = paste0('stage',5:20,'_')
  )
  sel_crit <- 'crit1'
}

# option: 2021alpha  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2021alpha'){
  ndays_sim <- 457
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('stage',20:29,'_'),
                                                 'log_VOC_alpha'),
                                    ignore_list = c('hosp',
                                                    'phi1_add')
  )
  sel_crit <- 'crit2'
}

# option: 2021delta  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2021delta'){
  ndays_sim <- 563 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('stage',29:39,'_'),
                                                 'log_VOC_delta'),
                                    ignore_list = c('hosp',
                                                    'phi1_add')
  )
  sel_crit <- 'crit2'
}

# option: 2021sept (september-november)  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2021sept'){
  ndays_sim <- 650 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c(paste0('stage',39:45,'_')),
                                    ignore_list = c('hosp',
                                                    'phi1_add')
  )
  sel_crit <- 'crit2'
}

# option: omicron (ba1ba2)----
if(exists('f_param') && !is.na(f_param) && grepl('omicron',f_param)){
  ndays_sim <- 730 
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_VOC_omicron_init',
                                                 'log_VOC_omicron_transm',
                                                 'log_VOC_omicron_gamma_factor',
                                                 paste0('stage',45:51,'_')))
  sel_crit <- 'crit2'
}

# option: hload (hosp and ICU load) ----
if(exists('f_param') && !is.na(f_param) && grepl('hload',f_param)){
  ndays_sim <- 730
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_delta3',
                                                 'phi1_add'),
                                    ignore_list = c('ba4ba5')

  )
  sel_crit <- 'crit3'
}

# option: mort (mortality) ----
if(exists('f_param') && !is.na(f_param) && grepl('mort',f_param)){
  ndays_sim <- 730
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('log_mu'),
                                    ignore_list = c('ba4ba5')
                                    
  )
  sel_crit <- 'crit4'
}

# option: 2020q1  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2020q1'){
  ndays_sim <- 183 
  
  opt_aggr_q_stages <- list(3:10)
  # account for aggregated q-param
  parms       <- adjust_q_param(parms=parms,parms_names=parms_names,
                                parms_names_estim = NA,
                                bool_setup = TRUE,
                                opt_waves = opt_aggr_q_stages)
  parms_names <- names(parms)

  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    sel_id = 1:min(which(grepl('log_n0',parms_names))),
                                    tag_list = c(paste0('stage_aggr',1,'_'),
                                                 paste0('stage',1:2,'_')),
                                    ignore_list = c('log_phi1',
                                                    'log_beta1',
                                                    'log_delta3',
                                                    'log_delta4',
                                                    'log_mu',
                                                    'mcmc')
  )
  length(parms_names_estim)
  sel_crit <- 'crit1'
}

# option: 2020q2  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2020q2'){
  ndays_sim <- 306 
  
  opt_aggr_q_stages <- list(3:10,17:19)
  # account for aggregated q-param
  parms       <- adjust_q_param(parms=parms,parms_names=parms_names,
                                parms_names_estim = NA,
                                bool_setup = TRUE,
                                opt_waves = opt_aggr_q_stages)
  parms_names <- names(parms)
  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    sel_id = 1:min(which(grepl('log_n0',parms_names))),
                                    tag_list = c(paste0('stage_aggr',1:2,'_'),
                                                 paste0('mat_w',c(1:2,11:16),'_')),
                                    ignore_list = c('log_phi1',
                                                    'log_beta1',
                                                    'log_delta3',
                                                    'log_delta4',
                                                    'log_mu',
                                                    'mcmc')
  )
  length(parms_names_estim)
  sel_crit <- 'crit1'
}

# option: 2020q3  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2020q3'){
  ndays_sim <- 457 
  
  opt_aggr_q_stages <- list(3:10,17:19,20:28)
  # account for aggregated q-param
  parms       <- adjust_q_param(parms=parms,parms_names=parms_names,
                                parms_names_estim = NA,
                                bool_setup = TRUE,
                                opt_waves = opt_aggr_q_stages)
  parms_names <- names(parms)
  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    sel_id = 1:min(which(grepl('log_n0',parms_names))),
                                    tag_list = c(paste0('stage_aggr',1:3,'_'),
                                                 'log_VOC_alpha',
                                                 paste0('mat_w',c(1:2,11:16),'_')),
                                    ignore_list = c('log_phi1',
                                                    'log_beta1',
                                                    'log_delta3',
                                                    'log_delta4',
                                                    'log_mu',
                                                    'phi1_add',
                                                    'hosp',
                                                    'mcmc')
  )
  length(parms_names_estim)
  sel_crit <- 'crit2'
}

# option: 2020q4  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2020q4'){
  ndays_sim <- 563 
  
  opt_aggr_q_stages <- list(3:10,17:19,20:38)
  # account for aggregated q-param
  parms       <- adjust_q_param(parms=parms,parms_names=parms_names,
                                parms_names_estim = NA,
                                bool_setup = TRUE,
                                opt_waves = opt_aggr_q_stages)
  parms_names <- names(parms)
  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    sel_id = 1:min(which(grepl('log_comix',parms_names))-1),
                                    tag_list = c(paste0('stage_aggr',1:3,'_'),
                                                 'log_VOC_alpha',
                                                 'log_VOC_delta',
                                                 paste0('mat_w',c(1:2,11:16),'_')),
                                    ignore_list = c('log_phi1',
                                                    'log_beta1',
                                                    'log_delta3',
                                                    'log_delta4',
                                                    'log_mu',
                                                    'phi1_add',
                                                    'hosp',
                                                    'mcmc')
  )
  length(parms_names_estim)
  sel_crit <- 'crit2'
}


# option: 2020q7  ----
if(exists('f_param') && !is.na(f_param) && f_param == '2020q7'){
  ndays_sim <- 730 
  
  opt_aggr_q_stages <- list(1:8,9:11,12:17,18:31,32:42)
  # account for aggregated q-param
  parms       <- adjust_q_param(parms=parms,parms_names=parms_names,
                                parms_names_estim = NA,
                                bool_setup = TRUE,
                                opt_waves = opt_aggr_q_stages)
  parms_names <- names(parms)
  
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    sel_id = 1:min(which(grepl('log_comix',parms_names))-1),
                                    tag_list = c(paste0('stage_aggr',1:length(opt_aggr_q_stages),'_'),
                                                 'log_VOC',
                                                 paste0('mat_w',c(1:2,11:16),'_')),
                                    ignore_list = c('log_phi1',
                                                    'log_beta1',
                                                    'log_delta3',
                                                    'log_delta4',
                                                    'log_mu',
                                                    'phi1_add',
                                                    'hosp',
                                                    'ba4ba5',
                                                    'mcmc')
  )
  length(parms_names_estim)
  sel_crit <- 'crit2'
}

# option: load_for (hospital load and admissions for covid)
if(exists('f_param') && !is.na(f_param) && f_param == 'load_for'){
  
  parms_names <- names(parms)
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    tag_list = c('for_covid_coef',
                                                'log_delta3',
                                                 'phi1_add',
                                                 'log_mu')
  )
  sel_crit <- 'crit6'
}


# option debug (using dummy set: default) ----
if(exists('f_param') && !is.na(f_param) && f_param == 'debug'){
  parms_names_estim <- get_colnames(parms_names = parms_names,
                                    sel_id = c(6,9,20,125,200),
                                    tag_list = c('log_q_stage1_age5')
  )
  sel_crit <- 'crit2'
  ndays_sim <- 100
}

# check parms_names_estim ----
# if f_param is not recognized and parms_names_estim did not change, end MCMC calibration
if(length(parms_names_estim) == num_param_all){
  warning("PARAMETER SELECTION TAG UNKNOWN")
  return(-1)
}

# set rng seed
tag_numeric  <- as.numeric(charToRaw(run_tag));
seed_numeric <- sum(tag_numeric*10*seq(length(tag_numeric)))
set.seed(sum(tag_numeric*1:length(tag_numeric)))

# update calibration model horizon (could be changed during parameter selection)
parms['ndays_calibration'] <- ndays_sim

# adjust output directory path with current time and MCMC options
output_dir <- init_output_dir(output_tag,paste0(format(Sys.time(),'%Y%m%d_%H%M%S'),'_d',ndays_sim,'_e',length(parms_names_estim),'_i',n_iter,'_n',n_repeat,'_p',n_period,'_',sel_crit,'_',run_tag))
print(output_dir)

# prepare vaccine uptake
# this is done once, to save time each iteration
vacc_files_all <- dir(vaccine_uptake_dir,pattern = paste0(f_region,'.*.csv'),full.names = T)
i_vac_file     <- vacc_files_all[grepl(vaccine_uptake_pattern,vacc_files_all)]
vacc_schedule  <- read.csv(i_vac_file, header = T)

# vaccine protection: 1st and 2nd dose
V_mat         <- get_vaccine_protection_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_A_', parms = parms )
V_mat_dose2   <- get_vaccine_protection_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_B_', parms = parms )
V_mat_booster <- get_vaccine_protection_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_E_', parms = parms )
V_mat_2ndbooster <- get_vaccine_protection_matrix(vaccine_uptake = vacc_schedule, dose_tag = '_F_', parms = parms )


# save parameters names
write.table(parms_names_estim,file= file.path(output_dir,'MCMC_parameter_estim.csv'),sep=',',row.names = F, col.names = F)
write.table(parms_names,file= file.path(output_dir,'MCMC_parameter_all.csv'),sep=',',row.names = F, col.names = F)

# set text file to log model output and MCMC progress
log_file_name <- file.path(output_dir,'MCMC_logfile.txt')
cat('\noutput tag:', file = log_file_name,append = T,fill = T)
cat(output_dir, file = log_file_name,append = T,fill = T)
cat('\nparameter priors:', file = log_file_name,append = T,fill = T)
cat(chains_param_file, file = log_file_name,append = T,fill = T)
cat('\nvaccine uptake:', file = log_file_name,append = T,fill = T)
cat(i_vac_file, file = log_file_name,append = T,fill = T)
cat(c('\nprevious paramer set:',f_chain_id), file = log_file_name,append = T,fill = T)
cat(c('\n',as.character(Sys.time())), file = log_file_name,append = T,fill = T)
cat(c('\nnumber of parameters:',length(parms_names_estim)), file = log_file_name,append = T,fill = T)

if(has_aggregated_q_param(parms_names_estim)){
  for(i_list in 1:length(opt_aggr_q_stages)){
    cat(c('\naggregated q-waves:',paste(opt_aggr_q_stages[[i_list]])), file = log_file_name,append = T,fill = T)
  }
}
cat('\n', file = log_file_name,append = T)
## Important functions
##-------------------- -
### Model 1: MCMC approach
###----------------------- -
log_likelihood <- function(parms, plots_code = "TRUE"){
  
    fit = log_likelihood_model(parms, 
                            contact_data = contact_data, 
                            method = "mean",
                            plots = plots_code,
                            parms_names = parms_names, 
                            ndays_sim = ndays_sim,
                            V_mat = V_mat,
                            V_mat_dose2 = V_mat_dose2,
                            V_mat_booster = V_mat_booster,
                            V_mat_2ndbooster = V_mat_2ndbooster,
                            be_ref_data = be_ref_data)

  return(list(crit                 = fit[[sel_crit]], 
              dev                  = -2*fit[[sel_crit]],
              new_hosp_since_intro = fit$new_hosp_since_intro,
              mean_new_hosp_icu    = fit$mean_new_hosp_icu,
              mean_new_deaths      = fit$mean_new_deaths,
              mean_new_discharged  = fit$mean_new_discharged))
}

log_prior <- function(parms){ 
  log_prior_model(parms,parms_names)
}


## Starting values ----
##---------------- -
startvalue        <- parms
startvalue[parms_names_estim] <- unlist(startvalue[parms_names_estim]) * runif(length(parms_names_estim),min = 1-s_deviance,max=1+s_deviance)

# check and limit adjusted q parameters
flag_M_coef <- grepl('log_q_stage',parms_names) 
startvalue[flag_M_coef & startvalue < -5] <- -4.9
startvalue[flag_M_coef & startvalue > 2]  <- 1.999

# check and limit adjusted transmission parameters
flag_transm <- grepl('log_VOC.*transm',parms_names)
startvalue[flag_transm & startvalue < 0] <- 0.01
startvalue[flag_transm & startvalue > 3.5]  <- 3.4

# check and limit adjusted gamma parameters
flag_gamma <- grepl('log_VOC.*gamma_factor',parms_names)
startvalue[flag_gamma & exp(startvalue) > 1e10] <- 1e9

# account for aggregated q-param
if(has_aggregated_q_param(parms_names_estim)){
  startvalue <- adjust_q_param(startvalue,parms_names,parms_names_estim,opt_waves = opt_aggr_q_stages)
}
# check and limit phi parameters
flag_phi <- grepl('log_phi',parms_names)
startvalue[flag_phi & startvalue < -13.8] <- -13
startvalue[flag_phi & startvalue > 4.59]  <- 4.5

# check and limit mu parameters
flag_mu <- grepl('log_mu',parms_names)
startvalue[flag_mu & startvalue < -13.8] <- -13
startvalue[flag_mu & startvalue > 4.59]  <- 4.5

# account for equal q parameters for different contact stages
startvalue <- adjust_q_param(startvalue,parms_names,parms_names_estim,opt_waves = opt_aggr_q_stages)

# startvalue['log_VOC_omicron_init'] <- 1
# startvalue['VOC_omicron_start'] <- 670

# check start values, and revert if new parameter value is not valid
scm_prior <- log_prior(startvalue)
if(length(scm_prior$invalid_param)>0){
  startvalue[scm_prior$invalid_param] <- parms[scm_prior$invalid_param]
  cat('\n\nreset invalid start values:',scm_prior$invalid_param, file = log_file_name,append = T)
}

# if still invalid values present from previous calibrations, use median of all prior values if the value is present
# note: median overcomes rounding issues
scm_prior <- log_prior(startvalue)
scm_prior$invalid_param <- scm_prior$invalid_param[scm_prior$invalid_param %in% names(parms_chains)]
if(length(scm_prior$invalid_param) == 1){
  startvalue[scm_prior$invalid_param] <- median(parms_chains[,scm_prior$invalid_param])
  cat('\n\nuse overall median for invalid start values:\n',scm_prior$invalid_param, file = log_file_name,append = T)
} else if(length(scm_prior$invalid_param) > 1){
  startvalue[scm_prior$invalid_param] <- apply(parms_chains[,scm_prior$invalid_param],2,median)
  cat('\n\nuse overall median for invalid start values:\n',scm_prior$invalid_param, file = log_file_name,append = T)
}

scm_prior2 <- log_prior(startvalue)
cat('\n\ninitial LL:\n',scm_prior2$log_prior_val, file = log_file_name,append = T)

# save start values
vec_chain_start     = file.path(output_dir,"MCMCstart_scm.csv")
write.table(t(startvalue), file = vec_chain_start, sep=",",row.names=F)

# save parameter info
write.table(data.frame(name= parms_names,
                       value = startvalue,
                       exp_value = exp(startvalue)),
            file= file.path(output_dir,'MCMC_parameter_values.csv'),sep=',',row.names = F, col.names = T)

## LaplacesDemon estimation (MCMC sampler)
##---------------------------------------- -
Model_CB = function(pars,Data){
  
  log_prior <- function(parms){ 
    log_prior_model(parms,parms_names)
  }
  ## select model parameters to calibrate
  pars_estim                    <- startvalue
  pars_estim[parms_names_estim] <- pars
  
  pars_estim <- adjust_q_param(pars_estim,parms_names,parms_names_estim,opt_waves = opt_aggr_q_stages)
  
  ## log(prior densities)
  LPr = log_prior(pars_estim)
  ## log-likelihood
  LL = log_likelihood(pars_estim, plots_code = "FALSE")
  ## log-posterior
  LP = LL$crit + LPr$log_prior_val
  ## Additional parameters
  R0 = LPr$R0
  
  ## output
  modelout = list(LP      = LP, 
                  Dev     = -2*LL$crit,
                  Monitor = c(pars,LP,R0,LPr$log_prior_val),
                  yhat    = LL$mean_new_hosp_icu, 
                  parm    = pars)
  return(modelout)
}

# compile R function to speed up MCMC deamon
library(compiler)
Model_CB   <- cmpfun(Model_CB)
J          <- length(parms_names_estim)       # Total number of model parameters + residual variance
mon.names  <- c(as.parm.names(list(beta=rep(0,J))),"LP","R0","LPr")
parm.names <- as.parm.names(list(beta=rep(0,J)))
MyData     <- list(J          = J, 
                   mon.names  = mon.names, 
                   parm.names = parm.names, 
                   y          = be_ref_data$hospital_admissions[1:ndays_sim])

# check initial parameter set
if(bool_plot_model_output){
  pdf(file.path(output_dir,'MCMC_SCM_start.pdf'),10,10)
  fit <- log_likelihood(as.numeric(startvalue),plots_code = TRUE)
  dev.off()
}

## Adaptive Metropolis-within-Gibbs algorithm (using the starting values)
##---------------------------------------------------------------------- -
ptm_sys <- Sys.time()
res_AMWG <- LaplacesDemon(Model          = Model_CB, 
                          Data           = MyData, 
                          Initial.Values = unlist(startvalue[parms_names_estim]),
                          Iterations     = n_iter, 
                          Status         = n_status, 
                          LogFile        = log_file_name,
                          Thinning       = 1,
                          Algorithm      = "AMWG", 
                          Specs          = list(B=NULL, n=n_repeat, Periodicity=n_period))

# add info to log file
cat(as.character(Sys.time()), file = log_file_name,append = T,fill = T)
cat('\n', file = log_file_name,append = T)

cat('MCMC Elapsed Time: ', file = log_file_name,append = T)
cat(as.character(format(Sys.time()-ptm_sys,digits=3)), file = log_file_name,append = T,fill = T)
cat('\n', file = log_file_name,append = T)


# Save Deamon
saveRDS(res_AMWG,file=file.path(output_dir,'MCMC_deamon.rds'))

## Save estimated parameters
##-------------------- -
chain_step1       = res_AMWG$Monitor[,1:length(parms_names_estim)]
add_chain_step1   = res_AMWG$Monitor[,length(parms_names_estim)+1]
dev_chain_step1   = res_AMWG$Monitor[,length(parms_names_estim)+2]
prior_chain_step1 = res_AMWG$Monitor[,length(parms_names_estim)+3]

## merge estimated and original parameters
chain_estim        <- data.frame(chain_step1)
names(chain_estim) <- parms_names_estim
chain_other        <- t(as.matrix(startvalue[!parms_names %in% parms_names_estim]))
chain_full         <- data.frame(chain_estim,chain_other)
dim(chain_full)

# option to account for aggregated q-param
chain_full <- adjust_q_param(chain_full,parms_names,parms_names,opt_waves = opt_aggr_q_stages)

# reorder columns (ordered by initial column names)
chain_full     <- chain_full[,parms_names]
chain_full_out <- as.matrix(chain_full)

## explore
if(bool_plot_model_output){
  pdf(file.path(output_dir,'MCMC_SCM_final_step.pdf'),10,10)
  fit <- log_likelihood(as.numeric(chain_full[nrow(chain_full),]),plots_code = TRUE)
  dev.off()
  
  plot(res_AMWG, 
       BurnIn=0, 
       MyData, 
       PDF=TRUE, 
       Parms=NULL, 
       FileName = file.path(output_dir,'MCMC_marginal_posterior_samples.pdf'))
  
  pdf(file.path(output_dir,'MCMC_caterpillar_plot.pdf'),10,10)
  caterpillar.plot(res_AMWG, Parms="beta")
  dev.off()
  
  cat('MCMC COMPLETE at ', file = log_file_name,append = T)
  cat(as.character(Sys.time()), file = log_file_name,append = T,fill = T)
  cat('\n', file = log_file_name,append = T)
}

# add MCMC info ----
chain_full$mcmc_chain_id       <- f_chain_id
chain_full$mcmc_iter_id        <- 1:nrow(chain_full)
chain_full$mcmc_ll_prior       <- prior_chain_step1
chain_full$mcmc_ll_posterior   <- add_chain_step1
chain_full$mcmc_core_version   <- get_scm_version()

# reorder columns to present MCMC info first
names_all    <- names(chain_full)
names_sorted <- c(names_all[grepl('mcmc_',names_all)],names_all[!grepl('mcmc_',names_all)])
chain_full   <- chain_full[,names_sorted]

# save chain as RDS
saveRDS(chain_full,file=file.path(output_dir,'MCMC_full.rds'))

# set filenames
vec_chain_step1       = file.path(output_dir,"MCMCchain_scm.csv")
vec_dev_chain_step1   = file.path(output_dir,"MCMCdevchain_scm.csv")

# store results
write.table(chain_full,file = vec_chain_step1, sep=",",row.names=F,col.names = T)
write.table(dev_chain_step1,file = vec_dev_chain_step1, sep=",",row.names=F,col.names = F)

# TMP: inspect this chain
merge_chains(output_dir=output_dir)

}

################################################################ #
## RUN THE MCMC FUNCTION ----
################################################################ #

# load this script, and run it
print("script loaded... start running")
if(length(args)==2){
  run(args[[1]],args[[2]])
} else if(length(args)==4){
  run(args[[1]],args[[2]],args[[3]],args[[4]])
} else if(length(args)==5){
  run(args[[1]],args[[2]],args[[3]],args[[4]],args[[5]])
} else if(length(args)==6){
  run(args[[1]],args[[2]],args[[3]],args[[4]],args[[5]],args[[6]])
} else if(length(args)==7){
  run(args[[1]],args[[2]],args[[3]],args[[4]],args[[5]],args[[6]],args[[7]])
} else {
  run()
}



