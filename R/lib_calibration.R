########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is needed to estimate model parameters using MCMC
#
# Copyright 2024, SIMID                                
########################################################################### #

# load functions and data ----
source('R/main_vaccination.R')
library(foreach)

merge_chains <- function(output_dir=NA){
  
  if(!exists('output_dir') || is.na(output_dir)){
    output_dir <- 'output/MCMC_SCM/debug'
  }
  
  # option to remove final '/'
  output_dir_new <- output_dir
  if("/" == substr(output_dir_new,nchar(output_dir_new),nchar(output_dir_new)+1)){
    output_dir_new <- paste0(output_dir_new,basename(output_dir_new))
  }
  
  output_dir_results <- paste0(output_dir_new,'_final')
  output_dir_pdf     <- paste0(output_dir_new,'_final_pdf')
  if(!dir.exists(output_dir_results)){
    dir.create(output_dir_results)
    dir.create(output_dir_pdf)
  }
  
  mcmc_files       <- dir(output_dir,pattern = 'MCMCchain_scm', recursive = T,full.names = T)
  parm_estim_files <- dir(output_dir,pattern = 'MCMC_parameter_estim', recursive = T,full.names = T)  
  parm_all_files   <- dir(output_dir,pattern = 'MCMC_parameter_all', recursive = T,full.names = T)  
  SCM_final_files  <- dir(output_dir,pattern = 'MCMC_SCM_final_step', recursive = T,full.names = T)
  log_files        <- dir(output_dir,pattern = 'MCMC_logfile', recursive = T,full.names = T)
  
  # generate mcmc_tag without hour stamp and chain id
  mcmc_tag_all <- unlist(strsplit(basename(dirname(mcmc_files[1])),split = '_'))[]
  mcmc_tag     <- paste(mcmc_tag_all[-c(2,length(mcmc_tag_all))],collapse='_')
  
  # initiate final tables
  mcmc_out      <- NULL
  mcmc_files <- data.frame(filename = mcmc_files,
                           mcmc_chain_id = NA)
  for(i_file in 1:nrow(mcmc_files)){
    s_file      <- mcmc_files$filename[i_file]
    chain_param <- read.table(s_file,header=T,sep=',')
    
    mcmc_files$mcmc_chain_id[i_file] <- unique(chain_param$mcmc_chain_id)
    mcmc_out                    <- rbind(mcmc_out,chain_param)
    
    print(c(mcmc_files$mcmc_chain_id[i_file],dim(chain_param)))
  }  
  
  # order mcmc_out by iteration and ll posterior
  mcmc_out   <- mcmc_out[order(mcmc_out$mcmc_iter_id,
                               mcmc_out$mcmc_ll_posterior,
                               decreasing = T),]
  
  # get parameter names
  dim(mcmc_out)
  param_names_estim <- unlist(read.table(parm_estim_files[1],header=F,sep=','))
  param_names       <- names(mcmc_out)
  
  # select n-last iterations
  n_iter           <- 5
  sel_iter         <- mcmc_out$mcmc_iter_id > (max(mcmc_out$mcmc_iter_id)-n_iter)
  
  # write to file
  chain_param <- write.table(mcmc_out[sel_iter,!is_aggregated_q_param(names(mcmc_out))], 
                             file=file.path(output_dir_results,
                                            paste0('MCMCmulti_',mcmc_tag,'.csv')),
                             sep=',',
                             row.names = F,
                             col.names = T)
  
  # write parameter summary ----
  tbl_final <- round(mcmc_out[mcmc_out$mcmc_iter_id == max(mcmc_out$mcmc_iter_id),],digits=3)
  tbl_start <- round(mcmc_out[mcmc_out$mcmc_iter_id == 1,],digits=3)
  chain_summary <- write.table(data.frame(parameter = param_names,
                                          singleton = as.numeric(colMeans(tbl_final) == apply(tbl_final,2,min)),
                                          selection = as.numeric(apply(tbl_start,2,min) != apply(tbl_final,2,min) | apply(tbl_start,2,max) != apply(tbl_final,2,max)),
                                          mean = colMeans(tbl_final),
                                          min = apply(tbl_final,2,min),
                                          max = apply(tbl_final,2,max),
                                          start_min = apply(tbl_start,2,min),
                                          start_max = apply(tbl_start,2,min)), 
                               file=file.path(output_dir_results,
                                              paste0('MCMC_parameter_values.csv')),
                               sep=',',
                               row.names = F,
                               col.names = T)
  
  # store summary file with posterior and prior log likelihood
  mcmc_ll_summary <- data.frame(mcmc_chain_id = tbl_final$mcmc_chain_id,
                                mcmc_ll_prior = round(tbl_final$mcmc_ll_prior),
                                mcmc_ll_posterior = round(tbl_final$mcmc_ll_posterior))
  mcmc_ll_summary  <- merge(mcmc_ll_summary,mcmc_files)
  mcmc_ll_summary  <- mcmc_ll_summary[order(mcmc_ll_summary$mcmc_ll_posterior),]
  chain_param <- write.table(mcmc_ll_summary,
                             file=file.path(output_dir_results,
                                            paste0('MCMC_parameter_summary.csv')),
                             sep=',',
                             row.names = F,
                             col.names = T)
  
  # copy final step files ----
  file.copy(SCM_final_files,
            file.path(output_dir_pdf,paste0(basename(dirname(SCM_final_files)),'.pdf')))
  
  # copy once param names + a log file ----
  file.copy(parm_estim_files[1],
            file.path(output_dir_results,paste0(basename(parm_estim_files[1]))))
  file.copy(parm_all_files[1],
            file.path(output_dir_results,paste0(basename(parm_all_files[1]))))
  file.copy(log_files[1],
            file.path(output_dir_results,paste0('chain_X_',basename(log_files[1]))))
  
  # explore estimated parameters
  mcmc_estim      <- mcmc_out[,param_names_estim]
  
  pdf(file=file.path(output_dir_results,paste0('chains_param_',mcmc_tag,'.pdf')))
  par(mfrow=c(3,3))
  i_param <- 1
  for(i_param in 1:length(param_names_estim)){
    plot(x = range(mcmc_out$mcmc_iter_id),
         y = range(pretty(range(mcmc_out[param_names_estim[i_param]]),2)),
         col=0,
         main = param_names_estim[i_param],
         ylab=param_names_estim[i_param],
         xlab='iteration')
    grid()
    for(i_chain in unique(mcmc_out$mcmc_chain_id)){
      lines(mcmc_out[mcmc_out$mcmc_chain_id == i_chain, param_names_estim[i_param]],
            col = i_chain+1,  # prevent color '0'
            lwd=2)
    }
  }
  dev.off()
  
  ## explore LL ----
  pdf(file=file.path(output_dir_results,paste0('chains_LL_',mcmc_tag,'.pdf')))
  
  # plot function
  plot_mcmc_ll_posterior  <- function(mcmc_out,y_select,plot_main=''){
    y_lim <- range(pretty(y_select))
    x_lim <- c(0,max(mcmc_out$mcmc_iter_id))
    plot(x_lim,
         y_lim,
         type='p',
         col=0,
         xlab='iteration',
         ylab='ll_posterior',
         main=plot_main)
    for(i_chain in unique(mcmc_out$mcmc_chain_id)){
      lines(mcmc_out$mcmc_iter_id[mcmc_out$mcmc_chain_id == i_chain],
            mcmc_out$mcmc_ll_posterior[mcmc_out$mcmc_chain_id == i_chain],
            col = i_chain+1)  # prevent color '0'
    }
  }
  
  mcmc_ll_posterior_final <- mcmc_out$mcmc_ll_posterior[mcmc_out$mcmc_iter_id == max(mcmc_out$mcmc_iter_id)]
  
  # final range
  plot_mcmc_ll_posterior(mcmc_out,mcmc_ll_posterior_final,'final range')
  
  # final top 50%
  plot_mcmc_ll_posterior(mcmc_out,quantile(mcmc_ll_posterior_final,c(0.5,1)),'final range')
  
  # overall top 20%
  plot_mcmc_ll_posterior(mcmc_out,quantile(mcmc_out$mcmc_ll_posterior,c(0.8,1)),'top 20%')
  
  # full range
  plot_mcmc_ll_posterior(mcmc_out,mcmc_out$mcmc_ll_posterior)
  
  dev.off()
}


