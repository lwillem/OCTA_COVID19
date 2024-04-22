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
output_tag <- 'wave1'

# get parameters: default
chains_param_files <- 'data/config/MCMCmulti_20230119_d730_belgium_ve30_customQ.csv'

# number of parameter sets, stochastic realisations and time horizon
num_param_sets      <- 5
num_stochastic_real <- 3
num_days_sim        <- 150

# set scenario tag
scen_tag    <- 'wave1'

################################################################ #
# MAIN: CALL MODEL ----
################################################################ #

# call function
projections_vaccination(output_tag=output_tag,
                        chains_param_files=chains_param_files,
                        scen_tag=scen_tag,
                        num_param_sets=num_param_sets,
                        num_stochastic_real=num_stochastic_real,
                        num_days_sim=num_days_sim
                        )


# print directory with model output
print(paste0("Inspect model output at './output/",output_tag,"'"))


