########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file can be used to run multiple simulation runs with different 
# social contact and/or vaccination uptake assumptions for Belgium or 
# one of the three regions. 
#
# Current parallelism: for different vaccine uptake file
#
# Copyright 2024, SIMID                                        
########################################################################### #

# clear workspace
rm(list=ls())

# load functions and data 
source('R/main_vaccination.R')

# set output tag (= directory name)
output_tag <- 'uptake'

# scenario tags and parameters (uncomment line to select one option)
scen_tag <- 'orig'       # to run with default VOC characteristics
# scen_tag <- 'omicron'   # to run without Omicron VOC

# get parameters: uncomment line to select VE_infectiousness (0 or 30%) and whether the q parameters are stage-specific ('customQ') or not ('aggrQ'). 
chains_param_files = dir('data/config',pattern='ve00_customQ.*.csv',full.names = T,recursive = T)  
# chains_param_files = dir('data/config',pattern='ve30_customQ.*.csv',full.names = T,recursive = T)
# chains_param_files = dir('data/config',pattern='ve00_aggrQ.*.csv',full.names = T,recursive = T) 
# chains_param_files = dir('data/config',pattern='ve30_aggrQ.*.csv',full.names = T,recursive = T) 

# number of parameter sets, stochastic realisations and time horizon
num_param_sets      <- 3
num_stochastic_real <- 2
num_days_sim        <- 730

#select vaccine uptake file(s) ----
vacc_files <- dir('data/uptake/uptake_vac_voc/',pattern = 'uptake_booster.*csv',full.names = T)

# option to adjust model parameters
adjusted_parameters <- data.frame()  # default
if(scen_tag == 'omicron'){
  adjusted_parameters <- data.frame(VOC_omicron_start = 1e3) 
}

# call function
projections_vaccination(output_tag          = output_tag,
                        chains_param_files  = chains_param_files,
                        scen_tag            = scen_tag,
                        num_param_sets      = num_param_sets,
                        num_stochastic_real = num_stochastic_real,
                        num_days_sim        = num_days_sim,
                        vacc_files          = vacc_files,
                        adjusted_parameters = adjusted_parameters)

