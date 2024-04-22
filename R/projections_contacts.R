################################################################################ #
# Stochastic Compartmental Model for SARS-CoV-2 in Belgium
# 
# This script is part of the Stochastic Compartmental Model project developed
# by the SIMID group. The project aims to model transmission dynamics of the
# SARS-CoV-2 virus in the context of Belgium during the COVID-19 pandemic.
#
# This particular file is designed to execute multiple simulation runs of the 
# COVID-19 spread within Belgium. The script can be executed as a standalone
# component
#
# Copyright 2024, SIMID                                     
################################################################################ #

# clear workspace
rm(list=ls())

# load packages, user-defined functions and data 
# !! make sure all packages from "R/install.R" are installed !!
source('R/main_vaccination.R')

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# set output tag (= directory name)
output_tag <- 'contacts'

# get parameters: default
chains_param_files <- 'data/config/MCMCmulti_20230120_d730_belgium_ve00_customQ.csv'

# number of model parameter sets, stochastic realisations and time horizon
num_param_sets      <- 1
num_stochastic_real <- 3
num_days_sim        <- 230

################################################################ #
# LOAD AND PRE-PROCESS ----
################################################################ #

# load reference data
be_ref_data <- get_latest_incidence_data()

# set vaccine uptake files (or NA if not applicable)
vacc_files <- NA

# set output directory
scen_tag    <- 'contacts'

# # load contact data
contact_data       <- readRDS('data/social_contact_data.rds')
contact_data$label <- 'CoMix'  # Note: This setting may be overwritten if scenario options 1 or 2 are activated
names(contact_data)

# Scenario option 1: re-use contact matrices
contact_data$db_C_sim <- rbind(contact_data$db_C_sim[1:6,],
                               contact_data$db_C_sim[1,])        # e.g. re-use contact data and q parameters form stage 1
contact_data$db_C_sim$date[nrow(contact_data$db_C_sim)] <- c('2020-06-15') # set start date of adjusted behaviour
contact_data$db_C_sim
contact_data$label <- 'prepandemic'  # Note: This setting may be overwritten if scenario options 2 is activated
summary(contact_data$db_C_sim)

# Scenario option 2: use new contact matrices
contact_data[['C_summer2020_asy']] <- contact_data$C_CoMix5_asy * 0.5    # include new matrix based on existing data but reduce contacts/transmission
contact_data[['C_summer2020_sy']]  <- contact_data$C_CoMix5_sy  * 0.25   # include new matrix based on existing data but reduce contacts/transmission
contact_data$db_C_sim              <- rbind(contact_data$db_C_sim,
                                            c('2020-07-01','C_summer2020',8))  # start on 1st of July with summer2020 matrix, and use q-parameters of stage 8
contact_data$label <- 'summer2020'
contact_data$db_C_sim
 
################################################################ #
# MAIN: CALL MODEL ----
################################################################ #

# call function
projections_vaccination(output_tag          = output_tag,
                        chains_param_files  = chains_param_files,
                        scen_tag            = scen_tag,
                        num_param_sets      = num_param_sets,
                        num_stochastic_real = num_stochastic_real,
                        num_days_sim        = num_days_sim,
                        be_ref_data         = be_ref_data,
                        contact_data        = contact_data
                        )

# print directory with model output
print(paste0("Inspect model output at './output/",output_tag,"'"))


################################################################ #
#' FYI: commands and help functions to run the modelling software
#' 
#'   sim_date2day('2020-05-01')        # help function: calendar date to model day index
#'   sim_day2date(61)                  # help function: model day index to calendar date
#'   get_M_change_day()                # contact behaviour change points (day index)
#'   sim_day2date(get_M_change_day())  # contact behaviour change points (calendar date)
#'   get_CoMix_change_day()            # CoMix contact behaviour change points (calendar date)
#'   get_CoMix_change_date()           # CoMix contact behaviour change points (calendar date)
#'   get_wave_colnames(wave_id = 1,bool_comix = TRUE)  # column names for first CoMix wave
#'   get_wave_colnames(wave_id = 1,bool_comix = FALSE) # column names for first additional wave

