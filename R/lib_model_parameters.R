########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to obtain, modify and explore model parameters.
#
# Copyright 2024, SIMID                    
########################################################################### #

update2latest_model_parameter_config <- function(parms){
  
  parms_names <- names(parms)
  
  if(length(parms_names)==0){
    stop("EMPTY PARAMETER LIST PROVIDED TO 'update2latest_model_parameter_config(params)'")
  }
  
  if(!any(grepl('log_VOC_delta_hr_hosp_age1',names(parms)))){
    names(parms) <- gsub('log_VOC_delta_hr_hosp','log_VOC_delta_hr_hosp_age1',names(parms))
    parms[paste0('log_VOC_delta_hr_hosp_age',2:10)] <- parms[grepl('log_VOC_delta_hr_hosp_age1',names(parms))]
    parms_names <- names(parms)
    print("include log_VOC_delta_hr_hosp_age",fill=T)
  }
  # 
  
  if(!any(grepl('p_asympt_age',names(parms)))){
    parms[paste0('p_asympt_age',1:10)] <- c(0.94,0.90,0.84,0.61,0.49,0.21,0.02,0.02,0.02,0.02)
    parms_names <- names(parms)
    print("include p_asympt_age",fill=T)
  }
  
  if(!any('f_asympt' %in% names(parms))){
    parms['f_asympt'] <- 0.51;   # relative infectiousness of asymptomatic vs. symptomatic cases
    parms_names <- names(parms)
    print("include f_asympt",fill=T)
  }
  
  if(!any(grepl('VOC_alpha',names(parms)))){
    sel_voc_alpha <- grepl('VOC_',parms_names) & !grepl('delta_',parms_names) & !grepl('omicron_',parms_names)
    parms_names[sel_voc_alpha] <- gsub('VOC_','VOC_alpha_',parms_names[sel_voc_alpha])
    names(parms) <- parms_names
    print("adjust parameter names to VOC_alpha",fill=T)
    
    sel_voc_delta3 <- grepl('delta3_',parms_names) & !grepl('wave2',parms_names)
    parms_names[sel_voc_delta3] <- gsub('delta3_','delta3_VOC_',parms_names[sel_voc_delta3])
    names(parms) <- parms_names
    print("adjust parameter names to delta3_VOC_",fill=T)
    
    sel_voc_mu <- grepl('mu_',parms_names) & !grepl('mu_age',parms_names) & !grepl('mu_sev',parms_names)
    parms_names[sel_voc_mu] <- gsub('voc','alpha',parms_names[sel_voc_mu])
    parms_names[sel_voc_mu] <- gsub('mu_','mu_VOC_',parms_names[sel_voc_mu])
    names(parms) <- parms_names
    print("adjust parameter names to mu_VOC_",fill=T)
    
    sel_voc_ve <- grepl('ve_',parms_names) & !grepl('delta',parms_names) & !grepl('omicron',parms_names) & !grepl('rate',parms_names) & !grepl('transmission',parms_names)
    parms_names[sel_voc_ve] <- gsub('ve_','ve_VOC_alpha_',parms_names[sel_voc_ve])
    parms_names <- gsub('ve_delta','ve_VOC_delta',parms_names)
    parms_names <- gsub('ve_omicron','ve_VOC_omicron',parms_names)
    names(parms) <- parms_names
    print("adjust parameter names to ve_VOC_",fill=T)
    
  }
  
  if(!any(grepl('VOC_ba4ba5',names(parms)))){
    sel_param_omicron <- parms[grepl('VOC_omicron',parms_names)]
    parms[gsub('omicron','ba4ba5',names(sel_param_omicron))] <- sel_param_omicron
    parms_names <- names(parms)
    
    parms['VOC_ba4ba5_init']  <- 0
    parms['VOC_ba4ba5_start'] <- 1e10 #760  # default: no introduction
    #parms['log_VOC_ba4ba5_transm'] <- parms['log_VOC_omicron_transm'] * 1.3
    print("include dummy parameter names for VOC_ba4ba5",fill=T)
  }
  
  if(!any(grepl('for_covid_coef',names(parms)))){
    parms['for_covid_coef'] <- 20000
    parms_names <- names(parms)
    print("default 'for_covid_coef parameter",fill=T)
  }
  
  if(!any(grepl('ve_.*_booster_waning',names(parms)))){
    #with waning booster similar to 2nd dose - infection only !!!  #temp
    parms['ve_waning_booster_rate'] <- 1/70 #parms['ve_waning_immunity_rate']
    parms['ve_VOC_alpha_infection_booster_waning']         <- parms['ve_VOC_alpha_infection_booster']        # reduction with respect to infection with waning immunity
    parms['ve_VOC_alpha_incr_severe_booster_waning'] <- parms['ve_VOC_alpha_incr_severe_booster']       # incremental reduction within dynamic model
    parms['ve_VOC_delta_infection_booster_waning']   <- parms['ve_VOC_delta_infection_booster']   # reduction with respect to infection with waning immunity
    parms['ve_VOC_delta_incr_severe_booster_waning'] <- parms['ve_VOC_delta_incr_severe_booster']
    parms['ve_VOC_omicron_infection_booster_waning']   <- parms['ve_VOC_omicron_infection_booster']   # reduction with respect to infection with waning immunity
    parms['ve_VOC_omicron_incr_severe_booster_waning'] <- parms['ve_VOC_omicron_incr_severe_booster']
    parms['ve_VOC_alpha_incr_severe_waning']       <- parms['ve_VOC_alpha_incr_severe_rna2']       # incremental reduction within dynamic model
    parms['ve_VOC_delta_incr_severe_waning'] <- parms['ve_VOC_delta_incr_severe_rna2']
    parms['ve_VOC_omicron_incr_severe_waning'] <- parms['ve_VOC_omicron_incr_severe_rna2']
    parms_names <- names(parms)
    print("custom ve_booster_waning parameters",fill=T)
  }
  
  if(!any(grepl('ve_.*_reinfvac',names(parms)))){
    parms['ve_waning_infection_rate'] <- 0 #1/730 #parms['ve_waning_booster_rate']
    parms['ve_VOC_alpha_infection_reinf']     <- parms['ve_VOC_alpha_infection_waning']
    parms['ve_VOC_alpha_incr_severe_reinf']   <- parms['ve_VOC_alpha_incr_severe_waning']
    parms['ve_VOC_delta_infection_reinf']     <- parms['ve_VOC_delta_infection_waning']
    parms['ve_VOC_delta_incr_severe_reinf']   <- parms['ve_VOC_delta_incr_severe_waning']
    parms['ve_VOC_omicron_infection_reinf']   <- parms['ve_VOC_omicron_infection_waning']
    parms['ve_VOC_omicron_incr_severe_reinf'] <- parms['ve_VOC_omicron_incr_severe_booster']
    
    parms['ve_waning_infection_booster_rate'] <- 0 # parms['ve_waning_booster_rate']
    parms['ve_VOC_alpha_infection_reinfvac']     <- parms['ve_VOC_alpha_infection_booster']
    parms['ve_VOC_alpha_incr_severe_reinfvac']   <- parms['ve_VOC_alpha_incr_severe_booster']   
    parms['ve_VOC_delta_infection_reinfvac']     <- parms['ve_VOC_delta_infection_booster']
    parms['ve_VOC_delta_incr_severe_reinfvac']   <- parms['ve_VOC_delta_incr_severe_booster']
    parms['ve_VOC_omicron_infection_reinfvac']   <- parms['ve_VOC_omicron_infection_booster']
    parms['ve_VOC_omicron_incr_severe_reinfvac'] <- parms['ve_VOC_omicron_incr_severe_booster']
    parms_names <- names(parms)
    print("default reinf parameters",fill=T)
  }
  
  if(!any(grepl('log_.*_hr_hosp',names(parms)))){
    parms['log_VOC_alpha_hr_hosp'] <- 1e-15
    parms['log_VOC_delta_hr_hosp'] <- log(2.26)
    parms[paste0('log_VOC_omicron_hr_hosp_age',1:10)] <- log(c(1.0,0.89,0.67,0.57,0.54,0.42,0.32,0.42,0.49,0.49))
    parms_names <- names(parms)
    print("default hr parameters",fill=T)
  }
  
  if(!any(grepl('log_spring2022_phi1_add',names(parms)))){
    parms['icu_spring2022_start']    <- sim_date2day('2022-03-20')
    parms['log_spring2022_phi1_add'] <- parms['log_VOC_omicron_phi1_add']
    parms_names <- names(parms)
    print("default spring22 parameters",fill=T)
  }
  
  if(!any(grepl('mu_VOC_alpha_sev',names(parms)))){
    colnames_mu2 <- names(parms)[grepl('mu2_sev',names(parms))]
    parms[gsub('mu2_sev','mu_VOC_alpha_sev',colnames_mu2)] <- parms[colnames_mu2]
    parms[gsub('mu2_sev','mu_VOC_delta_sev',colnames_mu2)] <- parms[colnames_mu2]
    parms_names <- names(parms)
    print("default VOC_alpha parameters",fill=T)
  }
  
  
  # convert Comix/Additional social contact data into "stages"
  if(!any(grepl('log_q_stage',parms_names))){
    q <- exp(parms['log_q'])
    
    CoMix_coef  = exp(parms[grepl('log_comix_coef',parms_names)])
    Add_coef    = exp(parms[grepl('log_add_coef',parms_names)])
    
    CoMix_coef_id <- as.numeric(gsub('_age.*','',gsub('log_comix_coef_w','',names(CoMix_coef))))
    Add_coef_id   <- as.numeric(gsub('_age.*','',gsub('log_add_coef_mat_w','',names(Add_coef))))
    
    q_by_stage <- c(rep(1,10),
                    Add_coef[Add_coef_id %in% 1:2],
                    CoMix_coef[CoMix_coef_id %in% 1:8],
                    Add_coef[Add_coef_id %in% 3:8],
                    CoMix_coef[CoMix_coef_id %in% 9:max(CoMix_coef_id)]
    ) * q
    
    log_q_by_stage <- log(q_by_stage)
    
    num_stages <- length(q_by_stage)/10
    names(log_q_by_stage) <- paste0('log_q_stage',rep(1:num_stages,each=10),'_age',rep(1:10,num_stages))
    length(log_q_by_stage)
    #names(log_q_by_stage)
    # 
    # include new "stage parameters" and remove previous notations
    parms <- parms[!grepl('log_.*_coef.*age',names(parms))]
    parms <- parms[!grepl('log_q',names(parms))]
    parms[names(log_q_by_stage)] <- log_q_by_stage
    
    print("make use of generic 'stage' parameters",fill=T)
  }

  if(!any(grepl('behavioral_change_delay',names(parms)))){
    parms['behavioral_change_delay'] <- 5
    print("make use of default 'behavioral_change_delay = 5'",fill=T)
  }

  return(parms)
}

aggregate_q_parameters <- function(parms_chains,chains_param_file=NA,chain_selection=NA){
  
  # option to provide parameter file name
  if(!is.na(chains_param_file)){
    if(missing(parms_chains)){
      parms_chains     = read.table(chains_param_file, sep = ",", header = T)
    } else{
      warning('The provided "chains_param_file" to "explore_q_param()" is overruled.')
    }
  } 
  
  if(!any(is.na(chain_selection))){
    parms_chains <- parms_chains[chain_selection,]
  }
  
  # set number of age groups and waves
  num_age_groups <- 10
  # nwave_total <- length(get_M_change_day())-1
  nwave_total <- identify_total_nb_stages(names(parms_chains))
  
  # i_param <- 1
  parm_names <- names(parms_chains)

  q_summary <- array(NA,dim=c(nwave_total,3,num_age_groups),
                         dimnames=list(list(),
                                       c('q_mean','q_LL','q_UL'),
                                       c(paste0('age',1:num_age_groups))))
  i_age <- 1
  for(i_age in 1:10){
    #parms_age <- exp(parms_chains[,parm_names[grepl(paste0('log_.*_coef_.*_age',i_age),parm_names)]])
    parms_age <- exp(parms_chains[,parm_names[grepl(paste0('log_q_stage.*_age',i_age),parm_names)]])
    parms_age <- parms_age[,gsub('log_q_stage.*_age','',names(parms_age)) == i_age]
    names(parms_age)
    dim(parms_age)
    
    q_summary[,,i_age] <- cbind(q_mean = apply(parms_age,2,mean),
                                q_LL = apply(parms_age,2,min),
                                q_UL = apply(parms_age,2,max))
  }
  
  return(q_summary)
}

explore_q_param <- function(parms_chains,db_C_sim,x_lim=NULL,chains_param_file=NA){

  num_age_groups <- 10

  # reduce db_C_sim to the given parameter set
  num_stages <- identify_total_nb_stages(names(parms_chains))
  sel_C_sim <- db_C_sim[1:num_stages,]
  sel_C_sim$date <- as.Date(sel_C_sim$date)
  
  q_summary <- aggregate_q_parameters(parms_chains=parms_chains)
  
  age_meta <- data.frame(tag = get_age_groups(),
                         # col = brewer.pal(10,'Set3'),
                         col = rainbow(num_age_groups),
                         #col = 1:10,
                         pch = 1:num_age_groups)
  
  ## Q (polygon)
  i_age <- 1
  par(mfrow=c(3,4))
  for(i_age in 1:10){
    y_max <- max(0.5,q_summary)
    plot(sel_C_sim$date,
         q_summary[,'q_mean',i_age],
         col=0,
         ylim = c(0,y_max),
         xaxt='n',
         main = paste('Age',age_meta$tag[i_age]),
         ylab='q',
         xlab='Date')
    grid(ny=NULL,nx=NA)
    add_polygon(scenario_data = q_summary[,,i_age],
                scenario_dates = sel_C_sim$date,
                col=2)
    add_date_axis(sel_C_sim$date[-1])
    
    # add stages
    stage_id_sel <- pretty(sel_C_sim$q_stage_id)
    abline(v=sel_C_sim$date[stage_id_sel],
           lty=2,
           col=8)
    text(stage_id_sel,
         x=sel_C_sim$date[stage_id_sel],
         y=y_max,
         col=1,
         cex = 0.75)
  }
}

identify_total_nb_stages <- function(parms_names){
  
  col_tag <- 'log_q_stage'
  
  nb_age_groups    <- sum(grepl(paste0(col_tag,1,'_'),parms_names))
  nb_columns       <- sum(grepl(col_tag,parms_names))
  nb_stages        <- nb_columns / nb_age_groups
  
  return(nb_stages)
}

adjust_stage_param <- function(f_parms, contact_data){
  
  nstages_q        <- identify_total_nb_stages(names(f_parms))
  nstages_contact  <- nrow(contact_data$db_C_sim)
  
  if(nstages_q > nstages_contact){
    f_parms <- select_q_param(f_parms,nstages_contact)
  } 
  
  while(nrow(contact_data$db_C_sim) > identify_total_nb_stages(names(f_parms))){
      f_parms <- include_new_stage(f_parms=f_parms)
  }
  
  return(f_parms)
}

select_q_param <- function(f_params,nstages){
  
  col_tag        <- 'log_q_stage'
  bool_include <- !grepl(col_tag,names(f_params))
  
  for(i_stage in 1:nstages){
    bool_include <- bool_include | grepl(paste0(col_tag,i_stage,'_'),names(f_params))
  }
  print(paste("WARNING: REDUCED Q PARAMETERS TO",nstages,"STAGE(S)"))
  
  return(f_params[bool_include])
}

include_new_stage <- function(f_parms){
  
  col_tag <- 'log_q_stage'
  
  nstages          <- identify_total_nb_stages(names(f_parms))
  num_param        <- length(names(f_parms))
  i_w_age10        <- which(grepl(paste0(col_tag,nstages,'_age10'),names(f_parms)))
  i_new            <- c(1:i_w_age10,         # param 1 till stage X-1
                        i_w_age10- (9:0))    # new stage
  
  if(i_w_age10 < num_param){
    i_new <- c(i_new,
               (i_w_age10+1):length(names(f_parms))) # final param
  }
  
  if(is.null(ncol(f_parms))){
    f_parms <- f_parms[i_new]
  } else{
    f_parms <- f_parms[,i_new]
  }
  
  i_col_prev <- i_w_age10- (9:0)
  i_col_new  <- i_w_age10+(1:10)
  names(f_parms)[i_col_new] <- gsub(paste0('stage',nstages,'_'),
                                    paste0('stage',nstages+1,'_'),
                                    names(f_parms[i_col_prev]))
  
  print(paste("WARNING: DUPLICATED Q PARAMETERS FOR STAGE",nstages+1))
  
  return(f_parms)
}

get_colnames <- function(parms_names,
                         tag_list = NA,
                         ignore_list=NULL,
                         sel_id = NA){
  
  bool_selection <- matrix(FALSE,length(parms_names))
  if(!any(is.na(sel_id))){
    bool_selection[sel_id] <- TRUE
  }
  if(!any(is.na(tag_list))){
    for(i_tag in tag_list){
      bool_selection[grepl(i_tag,parms_names)] <- TRUE
    }
  }
  if(!any(is.na(ignore_list))){
    for(i_tag in ignore_list){
      bool_selection[grepl(i_tag,parms_names)] <- FALSE
    }  
  }

  return(parms_names[bool_selection])
}

has_aggregated_q_param <- function(param_names){
  return(any(is_aggregated_q_param(param_names)))
}

is_aggregated_q_param <- function(param_names){
  return(grepl('stage_aggr',param_names))
}

get_ve_severe<- function(ve_infection,ve_incr_severe){
  return(1-(1-ve_infection)*(1-ve_incr_severe))
}
get_ve_incr_severe <- function(ve_infection,ve_severe){
  return(1-(1-ve_severe)/(1-ve_infection))
}

# ve_tag <- 've'
get_ve_infection <- function(parms,ve_tag){
  return(get_ve_param(parms,ve_type='infection',ve_tag))
}

get_ve_infection_matrix <- function(parms,ve_tag,num_age_groups){
  
  ve_infection <- get_ve_param(parms,ve_type='infection',ve_tag)
  ve_matrix <- matrix(rep(unlist(ve_infection),each=num_age_groups),ncol=length(ve_infection),byrow = F)
  colnames(ve_matrix) <- paste0('vac_',names(ve_infection))
  
  return(ve_matrix)
}


get_single_ve_param <- function(parms,ve_type,ve_tag,health_state){
  col_name <- paste(ve_tag,ve_type,health_state,sep='_')
  #print(col_name)
  if(col_name %in% names(parms)){
    return(parms[col_name])
  } else{
    return(col_name=NA)
  }
}

# generic function
get_ve_param <- function(parms,ve_type='infection',ve_tag){
  
  ve_param <- data.frame(adeno1 = get_single_ve_param(parms,ve_type,ve_tag,'adeno1'))
  ve_param$adeno2    <- get_single_ve_param(parms,ve_type,ve_tag,'adeno2')   # reduction with respect to infection after 2nd adeno-based dose
  ve_param$rna1      <- get_single_ve_param(parms,ve_type,ve_tag,'rna1')   # reduction with respect to infection after 1st mRNA dose 
  ve_param$rna2      <- get_single_ve_param(parms,ve_type,ve_tag,'rna2')   # reduction with respect to infection after 2nd mRNA dose 
  ve_param$waning    <- get_single_ve_param(parms,ve_type,ve_tag,'waning')   # reduction with respect to infection with waning immunity
  ve_param$booster   <- get_single_ve_param(parms,ve_type,ve_tag,'booster')   # reduction with respect to infection with booster dose
  ve_param$booster_waning    <- get_single_ve_param(parms,ve_type,ve_tag,'booster_waning')   # reduction with respect to infection with booster dose with waning immunity
  ve_param$reinf     <- get_single_ve_param(parms,ve_type,ve_tag,'reinf')   #reinfection path
  ve_param$reinfvac  <- get_single_ve_param(parms,ve_type,ve_tag,'reinfvac')   #reinfection path after vaccination
  
  return(ve_param)
}

explore_ve_parameters <- function(chains_param_file){

  # load parameter file
  cat(fill = T)
  cat(chains_param_file,fill=T)
  parms_chains     = read.table(chains_param_file, sep = ",", header = T)

  # assume all VE parameters constant
  parms <- parms_chains[1,]
  
  VOC_ve_infection <- data.frame(t(get_ve_infection(parms = parms,ve_tag='ve_VOC_alpha')),
                                 t(get_ve_infection(parms = parms,ve_tag='ve_VOC_delta')),
                                 t(get_ve_infection(parms = parms,ve_tag='ve_VOC_omicron')),
                                 t(get_ve_infection(parms = parms,ve_tag='ve_VOC_ba4ba5')))
  
  VOC_ve_severe <- data.frame(t(get_ve_param(parms = parms,ve_type= 'severe',ve_tag='ve_VOC_alpha')),
                                 t(get_ve_param(parms = parms,ve_type= 'severe',ve_tag='ve_VOC_delta')),
                                 t(get_ve_param(parms = parms,ve_type= 'severe',ve_tag='ve_VOC_omicron')),
                                 t(get_ve_param(parms = parms,ve_type= 'severe',ve_tag='ve_VOC_ba4ba5')))
  
  VOC_ve_incr_severe <- data.frame(t(get_ve_param(parms = parms,ve_type= 'incr_severe',ve_tag='ve_VOC_alpha')),
                              t(get_ve_param(parms = parms,ve_type= 'incr_severe',ve_tag='ve_VOC_delta')),
                              t(get_ve_param(parms = parms,ve_type= 'incr_severe',ve_tag='ve_VOC_omicron')),
                              t(get_ve_param(parms = parms,ve_type= 'incr_severe',ve_tag='ve_VOC_ba4ba5')))
  
  names(VOC_ve_incr_severe) <- names(VOC_ve_infection) <- names(VOC_ve_severe) <- c('VOC_alpha','VOC_delta','VOC_omicron','VOC_ba4ba5')
  row.names(VOC_ve_incr_severe)[1] <- row.names(VOC_ve_infection)[1] <- row.names(VOC_ve_severe)[1] <- "adeno1"
  cat("\nVE_infection",fill = T);print(round(VOC_ve_infection,digits=2))
  #cat("\nVE_severe",fill = T);print(VOC_ve_severe)
  cat("\nVE_incr_severe",fill = T);print(round(VOC_ve_incr_severe,digits=2))
  
  # check VE_severe using model input VE_infection and VE_severe
  effective_VOC_ve_severe <- get_ve_severe(ve_infection = VOC_ve_infection,
                                       ve_incr_severe = VOC_ve_incr_severe)
  
  cat("\nEFFECTIVE VE_severe",fill = T);print(round(effective_VOC_ve_severe,digits=2))
  
  waning_immunity_rates <- parms[grepl('waning_.*_rate',names(parms))]
  waning_immunity_period <- 1/waning_immunity_rates
  print(rbind(waning_immunity_period))
  
  # transmission
  cat("\nVE_transmission",fill = T);print(round(parms[grepl('ve_.*trans',names(parms))],digits=2))
  
  # reinfection
  #names(parms)[grepl('reinf',names(parms))]
  
  
}
  
save_vaccine_parameters <- function(){
  
  # VE references
  # Bernal 2021, NEJM : VE infection
  # Stowe2021 (pre-print): VE hospital admission
  # https://www.nature.com/articles/d41586-021-02054-z: VE transmission
  # Bernard 2021: pre-print on VOC omicron in England.
  
  vaccine_param <- data.frame(
    ve_VOC_alpha_infection_adeno1    = 0.49,  # reduction with respect to infection after 1st adeno-based dose
    ve_VOC_alpha_infection_adeno2    = 0.74,  # reduction with respect to infection after 2nd adeno-based dose
    ve_VOC_alpha_severe_adeno1       = 0.76,     # reduction with respect to severe infection with hospital admission after 1st adeno-based dose
    ve_VOC_alpha_severe_adeno2       = 0.86,     # reduction with respect to severe infection with hospital admission after 2nd adeno-based dose
    
    ve_VOC_alpha_infection_rna1      = 0.48,  # reduction with respect to infection after 1st mRNA dose #Bernal2021
    ve_VOC_alpha_infection_rna2      = 0.94,  # reduction with respect to infection after 2nd mRNA dose #Bernal2021
    ve_VOC_alpha_severe_rna1         = 0.83,     # reduction with respect to severe infection with hospital admission after 1st mRNA dose
    ve_VOC_alpha_severe_rna2         = 0.95,     # reduction with respect to severe infection with hospital admission after 2nd mRNA dose
    
    ve_VOC_alpha_infection_waning   = 0.627 ,  # [Andrews 2022] reduction with respect to infection after 2nd (mRNA) dose with waning immunity
    ve_VOC_alpha_severe_waning      = 0.917 ,  # [Andrews 2022] reduction with respect to severe infection with hospital admission 2nd (mRNA) dose with waning immunity
    
    ve_VOC_alpha_infection_booster   = 0.94 ,  # [copy of ve_VOC_alpha_infection_rna2] reduction with respect to infection after mRNA booster dose
    ve_VOC_alpha_severe_booster      = 0.95 ,  # [copy of ve_VOC_alpha_severe_rna2] reduction with respect to severe infection with hospital admission after mRNA booster dose

    ve_VOC_alpha_infection_booster_waning   = 0.889 , # [Andrews 2022] reduction with respect to infection after mRNA booster dose
    ve_VOC_alpha_severe_booster_waning      = max(0.95 *0.9,0.889),  # [copy of ve_VOC_alpha_severe_rna2] reduction with respect to severe infection with hospital admission after mRNA booster dose
    
    # protection delta variant
    ve_VOC_delta_infection_adeno1 = 0.429,       # 0.30 , # reduction with respect to infection after 1st adeno-based dose
    ve_VOC_delta_infection_adeno2 = 0.828,       # 0.67 ,  # reduction with respect to infection after 2nd adeno-based dose
    ve_VOC_delta_severe_adeno1    = 0.952 * 4/5, # 0.71 ,  # reduction with respect to severe infection with hospital admission after 1st adeno-based dose
    ve_VOC_delta_severe_adeno2    = 0.952,       # 0.92 ,  # reduction with respect to severe infection with hospital admission after 2nd adeno-based dose
    
    ve_VOC_delta_infection_rna1 = 0.723,         # 0.36 ,  # reduction with respect to infection after 1st mRNA dose
    ve_VOC_delta_infection_rna2 = 0.909,         # 0.88 ,  # reduction with respect to infection after 2nd mRNA dose
    ve_VOC_delta_severe_rna1    = 0.987 * 4/5,   # 0.94 ,  # reduction with respect to severe infection with hospital admission after 1st mRNA dose
    ve_VOC_delta_severe_rna2    = 0.987,          # 0.96 ,  # reduction with respect to severe infection with hospital admission after 2nd mRNA dose
    
    ve_VOC_delta_infection_waning   = 0.627 ,  # [Andrews 2022] reduction with respect to infection after 2nd (mRNA) dose with waning immunity
    ve_VOC_delta_severe_waning      = 0.917 ,  # [Andrews 2022] reduction with respect to severe infection with hospital admission 2nd (mRNA) dose with waning immunity
    
    ve_VOC_delta_infection_booster = 0.951, # 0.88,# 0.959 ,  # reduction with respect to infection after mRNA booster dose
    ve_VOC_delta_severe_booster    = 0.987, #0.996 ,  # reduction with respect to severe infection with hospital admission after mRNA booster dose
    
    ve_VOC_delta_infection_booster_waning   = 0.889 , # [Andrews 2022] reduction with respect to infection after mRNA booster dose
    ve_VOC_delta_severe_booster_waning      = max(0.987 *0.9,0.889) ,  # 90% of booster
    
    # protection against omicron variant 
    ve_VOC_omicron_infection_adeno1 = 0.177, # 0.129 , # reduction with respect to infection after 1st adeno-based dose
    ve_VOC_omicron_infection_adeno2 = 0.489, #0.190 ,  # reduction with respect to infection after 2nd adeno-based dose
    ve_VOC_omicron_severe_adeno1    = 0.81 * 4/5, # 0.497 ,  # reduction with respect to severe infection with hospital admission after 1st adeno-based dose
    ve_VOC_omicron_severe_adeno2    = 0.81,       # 0.600 ,  # reduction with respect to severe infection with hospital admission after 2nd adeno-based dose
    
    ve_VOC_omicron_infection_rna1 = 0.315,      # 0.187 ,  # reduction with respect to infection after 1st mRNA dose
    ve_VOC_omicron_infection_rna2 = 0.655,      # 0.241 ,  # reduction with respect to infection after 2nd mRNA dose
    ve_VOC_omicron_severe_rna1    = 0.81 * 4/5, # 0.596 ,  # reduction with respect to severe infection with hospital admission after 1st mRNA dose
    ve_VOC_omicron_severe_rna2    = 0.81,       # 0.668 ,  # reduction with respect to severe infection with hospital admission after 2nd mRNA dose
    
    ve_VOC_omicron_infection_waning   = 0.088 ,  # [Andrews 2022] reduction with respect to infection after 2nd (mRNA) dose with waning immunity
    ve_VOC_omicron_severe_waning      = 0.57 ,  # [CDC] reduction with respect to severe infection with hospital admission 2nd (mRNA) dose with waning immunity
    
    ve_VOC_omicron_infection_booster = 0.672 ,  # reduction with respect to infection after mRNA booster dose
    ve_VOC_omicron_severe_booster    = 0.90 ,  # reduction with respect to severe infection with hospital admission after mRNA booster dose
    
    ve_VOC_omicron_infection_booster_waning   = 0.457 , # [Andrews 2022] reduction with respect to infection after mRNA booster dose
    ve_VOC_omicron_severe_booster_waning      = 0.90 * 0.90 ,  # [copy of ve_severe_rna2] reduction with respect to severe infection with hospital admission after mRNA booster dose
    
    # protection against infectiousness/transmission
    ve_transmission = 0,# 0.45,
    
    # start protection
    delay_protection_rna1   = 21,
    delay_protection_rna2   = 7,
    delay_protection_adeno1 = 21,
    delay_protection_adeno2 = 7,
    
    ve_waning_immunity_rate      = 1/180,  # 6 months
    ve_waning_booster_rate       = 1/70    # 10 weeks
  )
  
  # reinfection (vac) 
  sel_param_2doses_waning <- vaccine_param[grepl('ve_VOC.*infection_waning',names(vaccine_param))]
  vaccine_param[gsub('waning','reinf',names(sel_param_2doses_waning))] <- sel_param_2doses_waning
  sel_param_booster <- vaccine_param[grepl('ve_VOC.*infection_booster',names(vaccine_param)) & !grepl('waning',names(vaccine_param))]
  vaccine_param[gsub('booster','reinfvac',names(sel_param_booster))] <- sel_param_booster
  vaccine_param$ve_waning_infection_rate <- 1/180
  
  # severe disease with reinfection
  sel_param_booster <- vaccine_param[grepl('ve_VOC.*severe_booster',names(vaccine_param)) & !grepl('waning',names(vaccine_param))]
  vaccine_param[gsub('booster','reinf',names(sel_param_booster))]   <- sel_param_booster * 0.8
  vaccine_param[gsub('booster','reinfvac',names(sel_param_booster))] <- sel_param_booster
  vaccine_param$ve_waning_infection_booster_rate <- 1/180
  
  # set incremental protection against severe disease  
  sel_ve_infection <- vaccine_param[grepl('ve_VOC_.*_infection',names(vaccine_param))]
  length(sel_ve_infection)
  for(sel_ve in names(sel_ve_infection)){
    vaccine_param[gsub('infection','incr_severe',sel_ve)]  <- get_ve_incr_severe(vaccine_param[sel_ve],vaccine_param[gsub('infection','severe',sel_ve)])
  }
  
  # make copy for BA4BA5
  sel_param_omicron <- vaccine_param[grepl('VOC_omicron',names(vaccine_param))]
  vaccine_param[gsub('omicron','ba4ba5',names(sel_param_omicron))] <- sel_param_omicron

  write.table(vaccine_param,
             file='output/vaccine_parameters_20220712.csv',
             sep=',',
             col.names=T,
             row.names=F)
  
  flag_infection <- grepl('infection',names(vaccine_param)) & !grepl('rate',names(vaccine_param))
  flag_severe    <- grepl('severe',names(vaccine_param)) & !grepl('incr',names(vaccine_param))
  flag_incr      <- grepl('incr',names(vaccine_param))
  
  flag_alpha   <- grepl('alpha',names(vaccine_param))
  flag_delta   <- grepl('delta',names(vaccine_param))
  flag_omicron <- grepl('omicron',names(vaccine_param))

  vaccine_param <- unlist(vaccine_param)  
  
  print(cbind(alpha_infection = vaccine_param[flag_infection & flag_alpha],
        alpha_severe = vaccine_param[flag_severe & flag_alpha],
        
        delta_infection = vaccine_param[flag_infection & flag_delta],
        delta_severe  = vaccine_param[flag_severe & flag_delta],
        
        omicron_infection = vaccine_param[flag_infection & flag_omicron],
        omicron_severe = vaccine_param[flag_severe & flag_omicron])
  )
  
}

explore_param <- function(chains_param_file){
  
  # read file  
  parms_chains     = read.table(chains_param_file, sep = ",", header = T)
  
  # set parameter names
  parms_names <- names(parms_chains)
  
  # plot all parameters
  pdf('output/param_all.pdf',9,9)
  par(mfrow=c(5,5))
  for(i_param in 1:ncol(parms_chains)){
    
    boxplot(parms_chains[,i_param],
            main = parms_names[i_param])
  }
  dev.off()
  
  # identify all parameters that are fixed
  parms_min <- apply(parms_chains,2,min)
  parms_max <- apply(parms_chains,2,max)
  print(which(parms_min == parms_max))
}


explore_disease_history_parameters <- function(chains_param_file,n_chains=NA,sel_chains = NA){
  
  # load parameter file
  cat(fill = T)
  cat(chains_param_file,fill=T)
  parms_chains     = read.table(chains_param_file, sep = ",", header = T)
  
  # make sure all latest parameter names and (placeholder) values are included
  parms_chains <- update2latest_model_parameter_config(parms_chains)
  
  if(all(!is.na(c(n_chains,sel_chains)))){
    warning('SET n_chains OR sel_chains')
    return(NULL)
  }
  
  if(!is.na(n_chains)){
    parms_chains <- parms_chains[1:n_chains,]
  }
  
  if(any(!is.na(sel_chains))){
    if(('mcmc_iter_id' %in% names(parms_chains))){
      parms_chains <- parms_chains[(parms_chains$mcmc_iter_id == max(parms_chains$mcmc_iter_id) &
                                      parms_chains$mcmc_chain_id %in% sel_chains),]
    } else{
      warning('MISSING: mcmc_iter_id')
      return(NULL)
    }

  }
  
  # get transition parameters
  gamma  = unlist(exp(parms_chains['log_gamma']));   # average length of the latency period - 2 days
  theta  = unlist(exp(parms_chains['log_theta']));   # infectious period at pre-symptomatic stage - 3.5 days
  delta1 = unlist(exp(parms_chains['log_delta1']));  # infectious period at asymptomatic stage - 3.5 days
  delta2 = unlist(exp(parms_chains['log_delta2']));  # infectious period at (mild) symptomatic stage - 3.5 days
  omega = colMeans(exp(parms_chains[grepl('log_omega',names(parms_chains))]));  # waiting time between symptom onset and hospitalization
  
  cat('exposed period:',mean(1/gamma),fill = T) # average length of the latency period
  cat('pre-sympt infectious period:',mean(1/theta),fill = T) # infectious period at pre-symptomatic stage
  cat('asympt infectious  period:',mean(1/delta1),fill = T) # infectious period at asymptomatic stage
  cat('mild sympt infectious period:',mean(1/delta2),fill = T) # infectious period at (mild) symptomatic stage
  cat('severe sympt infectious period before hosp admission:',(1/omega)) # waiting time between severe symptoms and hospitalization
  cat(fill = T)
  
  omicron_gamma_factor <- colMeans(exp(parms_chains['log_VOC_omicron_gamma_factor'])); 
  # omicron_transm       <- colMeans(exp(parms_chains['log_VOC_omicron_transm'])); 
  # delta_transm         <- colMeans(exp(parms_chains['log_VOC_delta_transm'])); 
  # alpha_transm         <- colMeans(exp(parms_chains['log_VOC_transm'])); 
  print(omicron_gamma_factor)
  
  # # check VOC part 1: increased transmission
  alpha_transm     <- 100*mean(exp(parms_chains[,'log_VOC_alpha_transm'])-1);       # additional transmissibility
  alpha_transm_CrI <-100*quantile(exp(parms_chains[,'log_VOC_alpha_transm'])-1,c(0.025,0.975));       # additional transmissibility
  
  # check VOC part 2: increased transmission
  delta_transm <- 100*mean(exp(parms_chains[,'log_VOC_delta_transm'])/exp(parms_chains[,'log_VOC_alpha_transm'])-1);       # additional transmissibility
  delta_transm_CrI <- 100*quantile(exp(parms_chains[,'log_VOC_delta_transm'])/exp(parms_chains[,'log_VOC_alpha_transm'])-1,c(0.025,0.975));       # additional transmissibility

  # check VOC part 3: increased transmission
  omicron_transm <- 100*mean(exp(parms_chains[,'log_VOC_omicron_transm'])/exp(parms_chains[,'log_VOC_delta_transm'])-1);       # additional transmissibility
  omicron_transm_CrI <- 100*quantile(exp(parms_chains[,'log_VOC_omicron_transm'])/exp(parms_chains[,'log_VOC_delta_transm'])-1,c(0.025,0.975));       # additional transmissibility
 
  # check VOC part 4: increased transmission
  if(!any(grepl('VOC_ba4ba5',names(parms_chains)))){
    parms_chains[,'log_VOC_ba4ba5_transm'] <- NA
    parms_chains[,'log_VOC_ba4ba5_init'] <- NA
    parms_chains[,'VOC_ba4ba5_start'] <- NA
    #parms['VOC_ba4ba5_start'] <- 1e10 #760  # default: no introduction
    #parms['log_VOC_ba4ba5_transm'] <- parms['log_VOC_omicron_transm'] * 1.3
  }
  
  ba4ba5_transm <- 100*mean(exp(parms_chains[,'log_VOC_ba4ba5_transm'])/exp(parms_chains[,'log_VOC_omicron_transm'])-1);       # additional transmissibility
  ba4ba5_transm_CrI <- 100*quantile(exp(parms_chains[,'log_VOC_ba4ba5_transm'])/exp(parms_chains[,'log_VOC_omicron_transm'])-1,c(0.025,0.975),na.rm=T);       # additional transmissibility
  
  omicron_start        <- colMeans(parms_chains['VOC_omicron_start'])
  omicron_start_range  <- range(parms_chains['VOC_omicron_start'])
  omicron_init         <- colMeans(exp(parms_chains['log_VOC_omicron_init']))
  omicron_init_range   <- range(exp(parms_chains['log_VOC_omicron_init']))
  
  ba4ba5_start        <- colMeans(parms_chains['VOC_ba4ba5_start'])
  ba4ba5_start_range  <- range(parms_chains['VOC_ba4ba5_start'])
  ba4ba5_init        <- colMeans(exp(parms_chains['log_VOC_ba4ba5_init']))
  ba4ba5_init_range  <- range(exp(parms_chains['log_VOC_ba4ba5_init']))
  
  cat(fill = T)
  cat('Alpha VOC transmission advantage (wrt Wuhan):\t',alpha_transm,'[',alpha_transm_CrI,']',fill = T)
  cat('Delta VOC transmission advantage (wrt Alpha):\t',delta_transm,'[',delta_transm_CrI,']',fill = T)
  cat('Omicron VOC transmission advantage (wrt Delta):\t',omicron_transm,'[',omicron_transm_CrI,']',fill = T)
  cat('BA4BA5 VOC transmission advantage (wrt Omicron):',ba4ba5_transm,'[',ba4ba5_transm_CrI,']',fill = T)
  cat('Omicron exposed period (days):\t',mean(1/(gamma*omicron_gamma_factor)),fill = T) # average length of the latency period
  cat(fill = T)
  
  cat('Omicron start:\t',format(sim_day2date(round(omicron_start))),'[',format(sim_day2date(round(omicron_start_range))),']',fill = T) # average length of the latency period
  cat('BA4BA5 start:\t',format(sim_day2date(round(ba4ba5_start))),'[',format(sim_day2date(round(ba4ba5_start_range))),']',fill = T) # average length of the latency period
  
  cat('Omicron init:\t',round(omicron_init),'[',round(omicron_init_range),']',fill = T) # average length of the latency period
  cat('BA4BA5 init:\t',round(ba4ba5_init),'[',round(ba4ba5_init_range),']',fill = T) # average length of the latency period
  
  
  # hospital hazard / probability
  phi0  = colMeans(expit(parms_chains[grepl('log_phi0',names(parms_chains))]));  # proportion of symptomatically infected with mild symptoms
  p_vec = colMeans(parms_chains[,grepl('p_asympt_age',names(parms_chains))]);    # proportion of asymptomatic cases
  h    = 0.25
  
  if(length(p_vec)==0){
    p_vec <- c(0.94,0.90,0.84,0.61,0.49,0.21,0.02,0.02,0.02,0.02)
  }
  
  parms <- unlist(parms_chains[1,])
  parms_names <- names(parms)
  delta2 = unlist(exp(parms['log_delta2']));  # infectious period at (mild) symptomatic stage - 3.5 days
  VOC_start  = parms[grepl('VOC_.*_start',parms_names)]
  VOC_name   = gsub('_start','',names(VOC_start))
  VOC_name  
  
  log_VOC_XXX_hosp <- NA
  sel_VOC              <- "VOC_delta"
  for(sel_VOC in VOC_name){
    if(sel_VOC == VOC_name[1]) {log_VOC_XXX_hosp <- NA}
    sel_VOC_param        <- parms[grepl(sel_VOC,parms_names)]
    names(sel_VOC_param) <- gsub(paste0(sel_VOC,'_'),'',names(sel_VOC_param))
    names(sel_VOC_param) <- gsub(paste0('_',sel_VOC),'',names(sel_VOC_param))
    #sel_VOC_param
    
    VOC_hr_hosp     <- exp(sel_VOC_param[grepl('log_hr_hosp',names(sel_VOC_param))]);   # hospital hazard ratio vs previous VOC
    log_VOC_XXX_hosp<- get_log_VOC_hosp(phi0,delta2,h, odds_ratio = VOC_hr_hosp, log_VOC_previous = log_VOC_XXX_hosp)
    VOC_hosp        <- unlist(expit(log_VOC_XXX_hosp))
    VOC_phi0        <- phi0 - (phi0 * VOC_hosp)
    VOC_phi0
  
   
    p_summary <- data.frame(rbind(round(1-p_vec,3),     # probability to be symptomatic
                                  round(1-phi0,3),      # probability Imild to Isevere
                                  round(VOC_hr_hosp,4),    # VOC adjustment factor
                                  round(1-VOC_phi0,3),  # probability Imild to Isevere 
                                  round((1-VOC_phi0)*(1-p_vec),3))) # probability Ia1 to Isevere
    rownames(p_summary) <- paste(sel_VOC,c('p_symptomatic','p_hospital','voc_hr','p_adj_voc','p_tot_hospital'))
    names(p_summary) <- gsub('.*_age','age',names(p_summary))
    print(p_summary)
    
  }
  
  # waning immunity assumptions
  sel_waning_param        <- parms[grepl("waning.*rate",parms_names)]
  names(sel_waning_param) <- paste0('1/',names(sel_waning_param))
  print(t(t(1/sel_waning_param)))
  
}

explore_invalid_param <- function(chains_param_file,n_chains=60){
  
  # load parameter file
  cat(fill = T)
  cat(chains_param_file,fill=T)
  parms_chains     = read.table(chains_param_file, sep = ",", header = T)
  dim(parms_chains)
  
  # loop over last n config, and check ll
  parms_names <- names(parms_chains)
  i <- nrow(parms_chains)
  for(i in 1:n_chains){
    invalid_param <- log_prior_model(unlist(parms_chains[i,]),parms_names)$invalid_param
    if(!is.null(invalid_param)){
      cat(c(i,invalid_param),fill = T)
    }  
  }
}

explore_param_table <- function(chains_param_file,n_chains=60){
  
  # load parameter file
  cat(fill = T)
  cat(chains_param_file,fill=T)
  parms_chains     = read.table(chains_param_file, sep = ",", header = T)
  dim(parms_chains)
  
  # loop over last n config, and check ll
  parms_names <- names(parms_chains)
  
  for(i in 1:n_chains){
    invalid_param <- log_prior_model(unlist(parms_chains[i,]),parms_names)$invalid_param
    if(!is.null(invalid_param)){
      cat(c(i,invalid_param),fill = T)
    }  
  }
}
