########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to load data and help functions for the stochastic model. 
# Variables are provided via the global R environment, which should be improved
# in future versions.
#
# Copyright 2024, SIMID                       
########################################################################### #

# Libraries
#---------- -

# CRAN packages
suppressPackageStartupMessages(library(LaplacesDemon))
suppressPackageStartupMessages(library(poisbinom))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(zoo))          # rollmean
suppressPackageStartupMessages(library(openxlsx)) 
suppressPackageStartupMessages(library(EpiEstim)) 
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))

# requires R (≥ 4.1.0)
if("EpiLPS" %in% installed.packages()){
  suppressPackageStartupMessages(library(EpiLPS)) 
}

# User defined package
suppressPackageStartupMessages(library(simid.rtools)) 

# load core and help functions
source('R/lib_model_core.R')
source('R/lib_combined_plot.R')
source('R/lib_model_parameters.R')
source('R/lib_stochastic_model.R')
source('R/lib_social_contacts.R')
source('R/lib_projections_vaccination.R')

## Helper functions  ----
##----------------- -
expit = function(eta){exp(eta)/(1+exp(eta))}
# logit = function(pp){log(pp/(1-pp))}

# defensive programming
logit = function(pp){
  pp <- ifelse(pp<1,pp,0.9999)
  return(log(pp/(1-pp)))
}

proposalfunction <- function(parms, sd = 0.005){
  return(rnorm(length(parms), mean = parms, sd = sd))
}

round_tot <- function(x, digits = 0){
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x,na.rm = T)) - sum(y,na.rm = T))
  y[indices] <- y[indices] + 1
  return(y / up)
}


## Variants of Concern (VOC)  ----
##------------------------ -

# set VOC types
global_lib_voc <- data.table(name = c('VOC_alpha',
                                      'VOC_delta',
                                      'VOC_omicron',
                                      'VOC_ba4ba5'),
                      start_date_ll = c("2020-11-06", "2021-04-05", "2021-11-01","2022-3-20"), # lower limit
                      start_date_ul = c("2020-12-26", "2021-07-14", "2021-12-06","2022-4-30")  # upper limit
                      )

global_lib_voc$model_strain <- rep(2:1,nrow(global_lib_voc))[1:nrow(global_lib_voc)]
global_lib_voc$start_ll <- sim_date2day(global_lib_voc$start_date_ll)
global_lib_voc$start_ul <- sim_date2day(global_lib_voc$start_date_ul)



