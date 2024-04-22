########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to install all CRAN packages and dependencies 
#
# Copyright 2024, SIMID                                       
########################################################################### #

# CRAN packages
install.packages('useful')        # compare.list
install.packages('EpiEstim')      # Rt calculations
install.packages('knitr')         # to convert markdown into pdf
install.packages('RColorBrewer')  # extra colors
install.packages('scales')        # to use transparent colors (cfr. credible intervals)
install.packages('data.table')    # to use the data.table format
install.packages('zoo')           # rollmean
install.packages('openxlsx')      # to read Excel files
install.packages('LaplacesDemon') # MCMC 
install.packages('poisbinom')     # poison binomial 
install.packages('socialmixr')    # social contact data analyses

# CRAN packages that require at least R version 4.1.0
if(as.numeric(R.Version()$major) >= 4 &&
   floor(as.numeric(R.Version()$minor)) >=1){
  install.packages('EpiLPS')
}

# to install the GitHub "simid.rtools" package for help functions on parallel computing etc.
install.packages("remotes") # to install packages from github 
remotes::install_github("lwillem/simid_rtools",force=F,quiet=T)

# to install the Belgian EQ5D package 
# requires R (≥ 4.1.0)
if(as.numeric(R.Version()$major) >= 4 &&
   floor(as.numeric(R.Version()$minor)) >=1){
  devtools::install_github("brechtdv/EQ5D.be")
}