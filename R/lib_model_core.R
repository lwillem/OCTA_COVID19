########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains the core function of the stochastic model, in combination
# with the function to calculate the log likelihood and log prior distributions. 
# This file also contains a function to plot rudimentary model output.
#
# Copyright 2024, SIMID                                        
########################################################################### #

if(0==1){ # debug parameters
  if(exists('run_param')) parms <- run_param 
  parms_names <- names(parms)
  ndays_sim <- 740;
  vaccine_uptake <- vacc_schedule # vacc_schedule_uptake
  contact_data = contact_data; plots = "FALSE"
  method = "mean"; # method = "stochastic"; 
  db_cnt_adjust = NULL
}
run_model <- function(parms, contact_data, method = "mean",
                   parms_names = NA, ndays_sim, 
                   db_cnt_adjust = NULL,
                   bool_full_output=TRUE, 
                   vaccine_uptake = NA,
                   be_ref_data,
                   V_mat, V_mat_dose2,
                   V_mat_booster,V_mat_2ndbooster){
  
  # (re-)store parameter names if they are not given, or the parameter matrix don't contain them in e.g. the MCMC deamon
  if(any(is.na(parms_names))){
    parms_names <- names(parms)
  } else {
    names(parms) <- parms_names
  }
  
  num_age_groups <- 10
  
  ## make sure that the model parameters are in line with latest requirements
  # as such, update colnames and add default values if needed
  ##-------------------------------------------------------- -
  parms <- update2latest_model_parameter_config(parms)
  parms_names <- names(parms)
  
  
  # ## General model parameters -----
  # ##-------------------------------------------------------- -
  delta2 = exp(parms['log_delta2']);  # infectious period at (mild) symptomatic stage - 3.5 days
  phi0  = expit(parms[grepl('log_phi0',parms_names)]);  # proportion of symptomatically infected with mild symptoms
  p_vec = parms[grepl('p_asympt_age',parms_names)];     # proportion of asymptomatic cases 
  omega = exp(parms[grepl('log_omega',parms_names)]);          # waiting time between symptom onset and hospitalization
  
  # update mortality from September 2020
  mortality_changepoint = parms[grepl('mortality_changepoint',parms_names)]
  mu2_param = parms[grepl('mu2_sev',parms_names)]
  mu2_sev   = c(0,expit(mu2_param[c(1,1:8)])); # age-specific hospital fatality ratio - probability of dying
  delta3_wave2 = exp(parms['log_delta3_wave2']);    # recovery rate from hospital (from sept 2021)
  
  # updated ICU/hosp proportion from fall 2021 and spring 2022
  phi1_fall2021      = expit(parms[grepl('log_phi1_age',parms_names)] * parms['log_fall2021_phi1_add']); # extra: proportion of severly infected with regular hospitalization
  icu_fall2021_start = parms['icu_fall2021_start']
  phi1_spring2022      = expit(parms[grepl('log_phi1_age',parms_names)] * parms['log_spring2022_phi1_add']); # extra: proportion of severly infected with regular hospitalization
  icu_spring2022_start = parms['icu_spring2022_start']
  
  ## Model parameters related to the social contact behaviour  -----
  ##---------------------------------------------------------------------------------- -
  nStage      <- sum(grepl('log_q_stage',parms_names)) / num_age_groups 
  q_stage     <- matrix(exp(parms[grepl('log_q_stage',parms_names)]), nrow = nStage, ncol = num_age_groups, byrow = T);
  
  
  # Delay in behavioural change
  # note: a delay of 5 means the original behaviour on t=0, starting with a transition 
  # from t+1, and full effect by t+5. For example: seq(0,1,length=6)
  behavioral_change_delay <- parms['behavioral_change_delay']
  
  ## Model initialization: population & region  ---- 
  ##------------------------- -
  # based on confirmed cases
  #TODO: make this independent of reference data
  aggr_dat    <- be_ref_data[,grepl('cases_',names(be_ref_data))]
  age_dist_cc <- apply(aggr_dat[1:13,], 2, sum)
  total_nr_cc <- sum(age_dist_cc)
  rel_freq_cc <- age_dist_cc/total_nr_cc
  
  n0 = exp(parms['log_n0'])
  imported_cases = round(rel_freq_cc*n0*(1/(1-p_vec)),0);     
  
  sel_region  <- get_region(parms['region_id'])
  cohort.size <- get_regional_pop(region=sel_region)
  
  ## Transmissibility parameters  ---- 
  ##-------------------------------- -
  f = parms['f_asympt'] ;   # 0.51,  relative infectiousness of asymptomatic vs. symptomatic cases
  
  ## Stochastic model parameters  ---- 
  ##-------------------------------- -
  h = parms['h'];             # resolution of the binomial chains
  
  ## Data augmentation  ---- 
  ##---------------------- -
  times = seq(0, ndays_sim-h, h)
  times_day <- floor(times)

  ## Initialization VOC ---- 
  ##----------------------- -
  VOC_start  = parms[grepl('VOC_.*_start',parms_names)]
  VOC_name   = gsub('_start','',names(VOC_start))
  VOC_name  
  
  if(!all(VOC_name %in% global_lib_voc$name)){
    stop("ERROR: UNKNOWN VOC NAME IN PARAMETER FILE")
  }
  
  if(any(duplicated(VOC_start))){
    VOC_start[duplicated(VOC_start)] <- VOC_start[duplicated(VOC_start)]+1
    warning("WARNING: duplicated VOC_start value, increased '",names(VOC_start)[duplicated(VOC_start)],"' by one day")
  }
  
  # VOC-specific transmission factor 
  transm_fctr        = 1                    # start with factor 1 (= 2020 wild type)
  transm_fctr_VOC    = 1                    # start with factor 1 (= 2020 wild type)
  
  # keep track of VOC-related adjustment of the hospital admission probability
  log_VOC_XXX_hosp <- matrix(NA,ncol=num_age_groups,nrow=length(VOC_name))
  for(i_VOC in 1:length(VOC_name)){

    sel_VOC_param        <- parms[grepl(VOC_name[i_VOC],parms_names)]
    log_VOC_previous     <- NA
    if(i_VOC > 1) { log_VOC_previous = log_VOC_XXX_hosp[i_VOC-1,] }
    
    VOC_hr_hosp              <- exp(sel_VOC_param[grepl('log_.*hr_hosp',names(sel_VOC_param))]);   # hospital hazard ratio vs previous VOC
    log_VOC_XXX_hosp[i_VOC,] <- get_log_VOC_hosp(phi0,delta2,h, odds_ratio = VOC_hr_hosp, log_VOC_previous = log_VOC_previous)
  }
  
  
  ## Regional simulation?  ---- 
  ##----------------------- -
  pop_be     <- get_regional_pop(region = 'belgium')
  pop_model  <- get_regional_pop(region = sel_region)
  pop_factor   <- sum(pop_model)/sum(pop_be)
  if(pop_factor <1){
    imported_cases  <- round(imported_cases * pop_factor)
  }  

  ## Vaccination parameters  ---- 
  ##--------------------------- -
  
  # get VE data in matrix format, to be aligned with the age-specific lambda
  # start with equal VE against infection for wild type and alpha variant
  ve_infection_matrix <- get_ve_infection_matrix(parms,ve_tag='ve_VOC_alpha',num_age_groups)
  ve_voc_infection_matrix <- ve_infection_matrix
  
  # protection against infectiousness/transmission
  ve_transmission <- parms['ve_transmission']  #temp put to zero
  
  # vaccine related waning immunity
  ve_waning_immunity_rate <- parms['ve_waning_immunity_rate'] * h
  ve_waning_booster_rate  <- parms['ve_waning_booster_rate'] * h
  
  # infection related waning immunity
  ve_waning_infection_rate <- parms['ve_waning_infection_rate'] * h
  ve_waning_infection_booster_rate <- parms['ve_waning_infection_booster_rate'] * h
  
  # check vaccine PROTECTION matrices: dose 1 and 2
  # this is the result of uptake AND delay until protection.
  if(any(!is.na(vaccine_uptake))){
    V_mat         <- get_vaccine_protection_matrix(vaccine_uptake = vaccine_uptake, dose_tag = '_A_', parms = parms, ndays_sim = ndays_sim )
    V_mat_dose2   <- get_vaccine_protection_matrix(vaccine_uptake = vaccine_uptake, dose_tag = '_B_', parms = parms, ndays_sim = ndays_sim )
    V_mat_booster <- get_vaccine_protection_matrix(vaccine_uptake = vaccine_uptake, dose_tag = '_E_', parms = parms, ndays_sim = ndays_sim )
    V_mat_2ndbooster <- get_vaccine_protection_matrix(vaccine_uptake = vaccine_uptake, dose_tag = '_F_', parms = parms, ndays_sim = ndays_sim )
  } else if(missing('V_mat') || missing('V_mat_dose2')){
    V_mat <- matrix(0,nrow=length(times),ncol=20)
    V_mat_dose2 <- matrix(0,nrow=length(times),ncol=20)
  }

  # check booster vaccine PROTECTION matrix
  if(missing('V_mat_booster')){
    V_mat_booster <- V_mat[,11:20]*0
    colnames(V_mat_booster) <- gsub('A','E',colnames(V_mat_booster))
  }
  # check 2nd booster vaccine PROTECTION matrix
  if(missing('V_mat_2ndbooster')){
    V_mat_2ndbooster <- V_mat[,11:20]*0
    colnames(V_mat_2ndbooster) <- gsub('A','F',colnames(V_mat_2ndbooster))
  }
  
  # specify time steps with new vaccine PROTECTION
  bool_n_uptake                             <- rowSums(V_mat) > 0
  bool_n_uptake[rowSums(V_mat_dose2) > 0]   <- TRUE
  bool_n_uptake[rowSums(V_mat_booster) > 0] <- TRUE
  bool_n_uptake[rowSums(V_mat_2ndbooster) > 0] <- TRUE
  if(length(bool_n_uptake)<length(times)){
    bool_n_uptake <- c(bool_n_uptake,rep(FALSE,length(times) - length(bool_n_uptake)))
  }

  ## Set population categories and array ----
  ##--------------------------------------------------- -
  pop_categories <- data.frame(name=c('orig',
                                      'vac_rna1','vac_rna2','vac_adeno1','vac_adeno2',
                                      'vac_waning','vac_booster','vac_booster_waning',
                                      'vac_reinf','vac_reinfvac'))
  pop_categories$vaccine_related <- grepl('vac',pop_categories$name)
  pop_categories$min_2doses      <- !(grepl('orig',pop_categories$name) | grepl('1',pop_categories$name))
  pop_categories[9,-1]           <- FALSE #TODO
  pop_categories$reinfection     <- grepl('reinf',pop_categories$name)
  
  
  # to store the current population 
  pop_matrix <- get_population_matrix(num_age_groups)
  pop_array  <- array(0,
                      dim=c(nrow(pop_matrix),ncol(pop_matrix),nrow(pop_categories)),
                      dimnames = list(NULL, colnames(pop_matrix),pop_categories$name))
  
  # define sub-populations that do matter for re-infection, without vaccination, and do not count the "original" twice 
  bool_non_reinf <- !grepl('reinf',pop_categories$name) & !grepl('orig',pop_categories$name)
  
  # initiate infections
  pop_array[,c('S','E'),'orig'] <- c(cohort.size - imported_cases, imported_cases);
  
  
  ## Transmission function (mean or stochastic)  ---- 
  ##--------------------------------------------------- -
  if (method == "mean"){   #floor added to to with stochastic one - temp removed
    get_new_transitions <- function(size_vec, new_prob_vec){
      return(size_vec*new_prob_vec)
    }
    bool_stochastic <- FALSE
  } else if (method == "stochastic"){  #floor added to avoid some NA  - temp removed
    get_new_transitions <- function(size_vec, new_prob_vec){
      if(any(size_vec%%1!=0)){browser()}
      return(rbinom(length(size_vec), size = size_vec, prob = new_prob_vec))
    }
    bool_stochastic <- TRUE
  } else{
    warning(paste("UNKNOWN method:",method))
  }
    
  ##################################################################################### #
  #  nfull_array: keep track of wild type and VOC compartments over time                      
  ##################################################################################### #
  nfull_array <- array(0,dim=c(ndays_sim,dim(pop_array)),
                       dimnames=list(NULL,
                                     NULL,
                                     dimnames(pop_array)[[2]],
                                     dimnames(pop_array)[[3]]));

  nfull_array[1,,,'orig'] <- pop_array[,,'orig']
  
  ## Transition rates  ----
  ##----------------------------- -
  # to store the population changes in one time step
  transition_col_names <- get_transition_names()
  transition_array <- array(0,
                            dim=c(num_age_groups,length(transition_col_names),nrow(pop_categories)),
                            dimnames = list(NULL, transition_col_names,pop_categories$name))
  
  # to store the probability/rates to change health state (for the stochastic and deterministic model, respectively)
  prob_array <- transition_array
  
  # adjust hospital admission probability to vaccine-related protection (based on Alpha VOC)
  # start with 2 identical copies for non-VOC and VOC (which is updated along the way)
  prob_matrix_vac      <- get_prob_array(parms,calendar_time = 0,log_VOC_XXX_hosp = log_VOC_XXX_hosp, pop_categories = pop_categories, ve_tag='ve_VOC_alpha')
  bool_voc_transitions <- c(rep(FALSE,dim(prob_matrix_vac)[2]),rep(TRUE,dim(prob_matrix_vac)[2]))
  prob_array[,bool_voc_transitions,]  <- prob_matrix_vac
  prob_array[,!bool_voc_transitions,] <- prob_matrix_vac

  # specify grouped column names for the transition and probability arrays
  col_origin    <- get_transition_origin()
  col_uptake    <- c('S','R','Rvoc')
  col_recovered <- c('R','Rvoc')
  col_no_voc    <- c(unique(col_origin[!grepl('voc',col_origin)]),'R','D')
  
  # Contact behaviour and transmission: beta   ----
  #----------------------------------------------- -
  beta_list_all <- calculate_all_beta_matrices(contact_matrices=contact_data, 
                                               q_stage,
                                               behavioral_change_delay,
                                               db_cnt_adjust)
  
  # adapt to regional simulation?
  if(pop_factor <1){
    sel_items <- which(grepl('sy',names(beta_list_all)))
    for(i in sel_items){
      beta_list_all[[i]] <- beta_list_all[[i]] * 1/pop_factor
    }
  }

  # run parameters
  bool_vaccine    <- FALSE
  bool_voc        <- FALSE

  # variables to log vaccine uptake (and related issues)
  testrealvacin <- testvacin <- testvacout <- rep(0,num_age_groups) 
  total_vac_rna1_remain    <- rep(0,num_age_groups)
  total_vac_rna2_remain    <- total_vac_rna1_remain
  total_vac_adeno1_remain  <- total_vac_rna1_remain
  total_vac_adeno2_remain  <- total_vac_rna1_remain
  total_vac_booster_remain <- total_vac_rna1_remain
  total_vac_2ndbooster_remain  <- total_vac_rna1_remain

  step_index <- 1
  # note: 
  # - each time step starts from the situation in "pop_array"
  # - the results of a time step are stored at nfull_array[step_index+1,....]
  ## LOOP OVER TIME STEPS ----
  for (step_index in 1:(length(times)-1)){
    #if(step_index == 2590) {cat('BREAK!'); break}
    #if(step_index == 4138) {cat('BREAK!'); browser() }
    
    calendar_time = times[step_index]
    #if(calendar_time == VOC_start[1]) {cat('BREAK!'); break }
    
    # update contact matrix?
    beta_out <- select_beta(calendar_time,
                           beta_list_all,
                           db_cnt_adjust)
    if(beta_out$bool_update){
      betas_asy <- beta_out$asy * f
      betas_sy  <- beta_out$sy
    }

    ## Update Mortality and Recovery      ----
    ##----------------------------- -
    if(calendar_time == mortality_changepoint){
      
      prob_array[,'t_D_icu',] =
        prob_array[,'t_D_hosp',] = 1 - exp(-h*delta3_wave2*mu2_sev)

      prob_array[,'t_R_icu',] =
        prob_array[,'t_R_hosp',] = 1 - exp(-h*delta3_wave2*(1-mu2_sev))
    }
   
     ## Updated ICU ratio for fall 2021       ----
    ##----------------------------- -
    if(calendar_time == icu_fall2021_start){
      # update severity, with respect to probability to be admitted to ICU
      prob_array[,'t_I_hosp',] <- 1 - exp(-h*  phi1_fall2021  *omega)
      prob_array[,'t_I_icu',] <- 1 - exp(-h*(1-phi1_fall2021)*omega)
    }
    
    ## Updated ICU ratio for spring 2022       ----
    ##----------------------------- -
    if(calendar_time == icu_spring2022_start){
      # update severity, with respect to probability to be admitted to ICU
      prob_array[,'t_Ivoc_hosp',]   <- 1 - exp(-h*  phi1_spring2022  *omega)
      prob_array[,'t_Ivoc_icu',]    <- 1 - exp(-h*(1-phi1_spring2022)*omega)
      
    }
    
    ## Update VOC         ----
    ##----------------------------- -
    if(any(calendar_time == VOC_start)){

      # update optimisation boolean
      bool_voc <- TRUE;
      
      ## 1. DEFINE VOC-SPECIFIC VALUES
      sel_VOC              <- VOC_name[VOC_start == calendar_time]
      sel_VOC_param        <- parms[grepl(sel_VOC,parms_names)]
      names(sel_VOC_param) <- gsub(paste0(sel_VOC,'_'),'',names(sel_VOC_param))
      names(sel_VOC_param) <- gsub(paste0('_',sel_VOC),'',names(sel_VOC_param))
      sel_VOC_param
      
      # VOC parameters  
      VOC_init        <- exp(sel_VOC_param['log_init']) * pop_factor # population factor for regional simulation
      VOC_init_age    <- round(VOC_init/num_age_groups, 0)
      
      VOC_prob_array <- get_prob_array(parms,calendar_time = calendar_time, log_VOC_XXX_hosp = log_VOC_XXX_hosp, pop_categories = pop_categories)
      
      ## 2. REGISTER VOC-SPECIFIC VALUES in STRAIN 1 OR 2
      if(global_lib_voc$model_strain[global_lib_voc$name == sel_VOC] == 1){
        # switch transmission potential
        transm_fctr <- exp(sel_VOC_param['log_transm']);   # additional VoC transmissibility
        
        # switch VE against infection for wild type strain
        ve_infection_matrix <- get_ve_infection_matrix(sel_VOC_param,ve_tag='ve',num_age_groups)
        
        # introduce cases
        pop_array[,'E','orig'] <- VOC_init_age
        
        # update transmission rates/probabilities
        prob_array[,!bool_voc_transitions,] <- VOC_prob_array
      
      } else {
        transm_fctr_VOC <- exp(sel_VOC_param['log_transm']);   # additional VoC transmissibility
        
        # switch VE against infection for voc strain
        ve_voc_infection_matrix <- get_ve_infection_matrix(sel_VOC_param,ve_tag='ve',num_age_groups)
        
        # introduce cases
        pop_array[,'Evoc','orig'] = VOC_init_age
        
        # update transmission rates/probabilities
        prob_array[,bool_voc_transitions,] <- VOC_prob_array
      }
    }

    ## Force of infection, with vaccine- and infection-related reduction in transmissibility  ---- 
    ##------------------------------------------------------------------------------------------------------- -
    if(bool_voc | bool_vaccine){  # temp : previously infect with reduced infectiousness or not ???
      tmp_aggr <- sum_by_pop_type(pop_array[,,!pop_categories$vaccine_related]) + 
                   sum_by_pop_type(pop_array[,,pop_categories$vaccine_related])*(1-ve_transmission)

      lambda_t = transm_fctr*(betas_asy %*% (tmp_aggr[,'I_presym'] + tmp_aggr[,'I_asym']) +
                                betas_sy  %*% (tmp_aggr[,'I_mild'] + tmp_aggr[,'I_sev']));

      lambda_t_voc = transm_fctr_VOC*(betas_asy %*% (tmp_aggr[,'Ivoc_presym'] + tmp_aggr[,'Ivoc_asym'] ) +
                                        betas_sy  %*% (tmp_aggr[,'Ivoc_mild'] + tmp_aggr[,'Ivoc_sev']));
      
     } else{
      tmp_aggr <- sum_by_pop_type(pop_array[,,!pop_categories$vaccine_related])
      
      lambda_t = betas_asy %*% (tmp_aggr[,'I_presym'] + tmp_aggr[,'I_asym']) +
                   betas_sy  %*% (tmp_aggr[,'I_mild'] + tmp_aggr[,'I_sev']);
      
      lambda_t_voc <- lambda_t*0
    }
    
    # perform a recurrent operation once
    neg_h_lambda_t     <- -h*lambda_t
    neg_h_lambda_t_voc <- -h*lambda_t_voc
    
    ## Chain binomial sequences  ---- 
    ##----------------------------- -
    # infections
    prob_array[,'t_E','orig']    <- 1 - exp(neg_h_lambda_t)
    prob_array[,'t_Evoc','orig'] <- 1 - exp(neg_h_lambda_t_voc)
    transition_array[,,'orig']   <- get_new_transitions(pop_array[,col_origin,'orig'],prob_array[,,'orig'])
    
    # reinfections
    prob_array[,'t_E','vac_reinf']    <- 1 - exp((1-ve_infection_matrix[,'vac_reinf'])*neg_h_lambda_t)
    prob_array[,'t_Evoc','vac_reinf'] <- 1 - exp((1-ve_voc_infection_matrix[,'vac_reinf'])*neg_h_lambda_t_voc)
    transition_array[,,'vac_reinf']   <- get_new_transitions(pop_array[,col_origin,'vac_reinf'],prob_array[,,'vac_reinf'])
   
    if(bool_vaccine){

      t_E    <- 1-exp((1-ve_infection_matrix)*neg_h_lambda_t[,rep(1,9)])
      t_Evoc <- 1-exp((1-ve_voc_infection_matrix)*neg_h_lambda_t_voc[,rep(1,9)])
      
      prob_array[,'t_E',colnames(ve_infection_matrix)]    <- t_E
      prob_array[,'t_Evoc',colnames(ve_infection_matrix)] <- t_Evoc
      
      for(i_pop_type in pop_categories$name[pop_categories$vaccine_related]){
        transition_array[,,i_pop_type] <- get_new_transitions(pop_array[,col_origin,i_pop_type],
                                                              prob_array[,,i_pop_type])
      }
    }  

                                     
    ### 1. Flows between disease states  ---- 
    ###------------------------------------ -
    for(i_pop_type in pop_categories$name[!pop_categories$vaccine_related]){
      pop_array[,,i_pop_type] <- update_compartments_matrix(pop_array[,,i_pop_type],
                                                            transition_array[,,i_pop_type])
    }
    
    if(bool_vaccine){
      for(i_pop_type in pop_categories$name[pop_categories$vaccine_related]){
        pop_array[,,i_pop_type] <- update_compartments_matrix(pop_array[,,i_pop_type],
                                                              transition_array[,,i_pop_type])
      }
    }
    
    pop_array[pop_array<0] <- 0
    
    ### 2. Vaccination events: uptake and waning immunity  ---- 
    ###-------------------------- -
   # always true once there is any vaccine uptake to account for waning immunity
   if(bool_vaccine | bool_n_uptake[step_index]){
      
     # set vaccine boolean to TRUE
      bool_vaccine <- TRUE
      
      coef_vac <- 1
      total_vac_adeno1 = V_mat[step_index, 1:10]*coef_vac;
      total_vac_rna1   = V_mat[step_index, 11:20]*coef_vac;
      total_vac_adeno2 = V_mat_dose2[step_index, 1:10]*coef_vac;
      total_vac_rna2   = V_mat_dose2[step_index, 11:20]*coef_vac;
      total_vac_booster = V_mat_booster[step_index, 1:10]*coef_vac;
      total_vac_2ndbooster = V_mat_2ndbooster[step_index, 1:10]*coef_vac;
      
      # create temporary reduced vectors
      tmp_c_vac            <- pop_array[,col_uptake,'orig']
      tmp_vac_adeno1_c_vac <- pop_array[,col_uptake,'vac_adeno1']
      tmp_vac_rna1_c_vac   <- pop_array[,col_uptake,'vac_rna1']
      tmp_vac_adeno2_c_vac <- pop_array[,col_uptake,'vac_adeno2']
      tmp_vac_rna2_c_vac   <- pop_array[,col_uptake,'vac_rna2']
      tmp_vac_waning_c_vac <- pop_array[,col_uptake,'vac_waning']
      tmp_vac_booster_c_vac<- pop_array[,col_uptake,'vac_booster']
      tmp_vac_booster_waning_c_vac<- pop_array[,col_uptake,'vac_booster_waning']
      tmp_vac_reinf_c_vac  <- pop_array[,col_uptake,'vac_reinf']
      tmp_vac_reinfvac_c_vac<- pop_array[,col_uptake,'vac_reinfvac']
      
      # # calculate vaccination factor
      # nvac must account for individuals both in tmp_c_vac and tmp_vac_reinf_c_vac
      factor_alive_nvac        = 1/(rowSums(tmp_c_vac)+rowSums(tmp_vac_reinf_c_vac))
      factor_alive_vac_adeno1  = 1/rowSums(tmp_vac_adeno1_c_vac,)
      factor_alive_vac_rna1    = 1/rowSums(tmp_vac_rna1_c_vac)
      factor_alive_vac_waning  = 1/rowSums(tmp_vac_waning_c_vac+tmp_vac_adeno2_c_vac+tmp_vac_rna2_c_vac)
      factor_alive_vac_for2ndbooster  = 1/rowSums(tmp_vac_booster_waning_c_vac+tmp_vac_reinfvac_c_vac)
      
      # calculate vaccination factor
      factor_alive_nvac[factor_alive_nvac==Inf] <- 0    
      factor_alive_vac_adeno1[factor_alive_vac_adeno1==Inf] <- 0
      factor_alive_vac_rna1[factor_alive_vac_rna1==Inf] <- 0
      factor_alive_vac_waning[factor_alive_vac_waning==Inf] <- 0
      factor_alive_vac_for2ndbooster[factor_alive_vac_for2ndbooster==Inf] <- 0
      
      #### mRNA-based uptake  ----  
      ####---------------------------------- -
      
      new_vac_rna1 <- floor(tmp_c_vac * ((total_vac_rna1+total_vac_rna1_remain) * factor_alive_nvac))
      new_vac_rna1_reinf  <- floor(tmp_vac_reinf_c_vac * ((total_vac_rna1+total_vac_rna1_remain) * factor_alive_nvac))
      total_vac_rna1_remain =  total_vac_rna1_remain + total_vac_rna1 - (rowSums(new_vac_rna1)+rowSums(new_vac_rna1_reinf))
      new_vac_rna2  <- (floor(tmp_vac_rna1_c_vac * ((total_vac_rna2+total_vac_rna2_remain) * factor_alive_vac_rna1)))
      total_vac_rna2_remain =  total_vac_rna2_remain + total_vac_rna2 - rowSums(new_vac_rna2)

      #### Adeno-based uptake   ----  
      ####----------------------------------- -
      new_vac_adeno1 <- floor(tmp_c_vac * ((total_vac_adeno1+total_vac_adeno1_remain) * factor_alive_nvac))
      new_vac_adeno1_reinf  <- floor(tmp_vac_reinf_c_vac * ((total_vac_adeno1+total_vac_adeno1_remain) * factor_alive_nvac))
      total_vac_adeno1_remain =  total_vac_adeno1_remain + total_vac_adeno1 - (rowSums(new_vac_adeno1)+rowSums(new_vac_adeno1_reinf))
      new_vac_adeno2  <- floor(tmp_vac_adeno1_c_vac * ((total_vac_adeno2+total_vac_adeno2_remain) * factor_alive_vac_adeno1))
      total_vac_adeno2_remain =  total_vac_adeno2_remain + total_vac_adeno2 - rowSums(new_vac_adeno2)
      
      #### booster uptake   ----  
      ####----------------------------------- -
      new_vac_booster_w    <- floor(tmp_vac_waning_c_vac * ((total_vac_booster+total_vac_booster_remain) * factor_alive_vac_waning))
      new_vac_booster_ad    <- floor(tmp_vac_adeno2_c_vac * ((total_vac_booster+total_vac_booster_remain) * factor_alive_vac_waning))
      new_vac_booster_rna    <- floor(tmp_vac_rna2_c_vac * ((total_vac_booster+total_vac_booster_remain) * factor_alive_vac_waning))
      total_vac_booster_remain =  total_vac_booster_remain + total_vac_booster - rowSums(matrix(new_vac_booster_w+new_vac_booster_ad+new_vac_booster_rna,nrow=num_age_groups))
 
      #### 2ndbooster uptake   ----  
      ####----------------------------------- -
      new_vac_2ndbooster_w    <- floor(tmp_vac_booster_waning_c_vac * ((total_vac_2ndbooster+total_vac_2ndbooster_remain) * factor_alive_vac_for2ndbooster))
      new_vac_2ndbooster_rinf    <- floor(tmp_vac_reinfvac_c_vac * ((total_vac_2ndbooster+total_vac_2ndbooster_remain) * factor_alive_vac_for2ndbooster))
      total_vac_2ndbooster_remain =  total_vac_2ndbooster_remain + total_vac_2ndbooster - rowSums(matrix(new_vac_2ndbooster_w+new_vac_2ndbooster_rinf,nrow=num_age_groups))
      
      
      #### Waning immunity and booster dose uptake  ----  using get_new_transitions functions.
      ####----------------------------------- -
      new_waning_rna2    <- get_new_transitions(tmp_vac_rna2_c_vac,ve_waning_immunity_rate)
      new_waning_adeno2  <- get_new_transitions(tmp_vac_adeno2_c_vac,ve_waning_immunity_rate)
      new_waning_booster <- get_new_transitions(tmp_vac_booster_c_vac,ve_waning_booster_rate)
      
    
      #### Update vaccine compartments
      pop_array[,col_uptake,'orig']       <- tmp_c_vac - new_vac_rna1 - new_vac_adeno1
      pop_array[,col_uptake,'vac_reinf']  <- tmp_vac_reinf_c_vac - new_vac_rna1_reinf - new_vac_adeno1_reinf
      pop_array[,col_uptake,'vac_rna1']   <- tmp_vac_rna1_c_vac + new_vac_rna1 + new_vac_rna1_reinf     - new_vac_rna2
      pop_array[,col_uptake,'vac_adeno1'] <- tmp_vac_adeno1_c_vac + c(new_vac_adeno1) + c(new_vac_adeno1_reinf) - new_vac_adeno2
      
      pop_array[,col_uptake,'vac_rna2']    <- tmp_vac_rna2_c_vac + new_vac_rna2 - new_vac_booster_rna - new_waning_rna2 
      pop_array[,col_uptake,'vac_adeno2']  <- tmp_vac_adeno2_c_vac + new_vac_adeno2 - new_vac_booster_ad - new_waning_adeno2

      pop_array[,col_uptake,'vac_waning']  <- tmp_vac_waning_c_vac + new_waning_rna2 + new_waning_adeno2 - new_vac_booster_w
      pop_array[,col_uptake,'vac_booster'] <- tmp_vac_booster_c_vac + new_vac_booster_w + new_vac_booster_ad + new_vac_booster_rna - new_waning_booster + new_vac_2ndbooster_w + new_vac_2ndbooster_rinf
      pop_array[,col_uptake,'vac_booster_waning'] <- tmp_vac_booster_waning_c_vac + new_waning_booster - new_vac_2ndbooster_w
      pop_array[,col_uptake,'vac_reinfvac'] <-  tmp_vac_reinfvac_c_vac - new_vac_2ndbooster_rinf
     

      #check negative doses due to vaccination, and correct for overflow if needed
      if(any(pop_array[,col_uptake,'orig']<0)){
        uptake_vector <- pop_array[,col_uptake,'orig']
        uptake_vector[uptake_vector>0] <- 0 # set all non-negative values to 0
        pop_array[,col_uptake,'vac_rna2'] <- pop_array[,col_uptake,'vac_rna2'] + uptake_vector # correct for the negative values
        pop_array[pop_array<0] <- 0
      }
      
      if(any(pop_array[,col_uptake,'vac_waning']<0)){
        uptake_vector <- pop_array[,col_uptake,'vac_waning']
        uptake_vector[uptake_vector>0] <- 0 # set all non-negative values to 0
        pop_array[,col_uptake,'vac_adeno2'] <- pop_array[,col_uptake,'vac_adeno2'] + uptake_vector   # correct for the negative values
        pop_array[,col_uptake,'vac_waning'] <- pop_array[,col_uptake,'vac_waning'] - uptake_vector # - negative values == 0
      }
      
      if(any(pop_array[,col_uptake,'vac_adeno2']<0)){
        uptake_vector <- pop_array[,col_uptake,'vac_adeno2']
        uptake_vector[uptake_vector>0] <- 0 # set all non-negative values to 0
        pop_array[,col_uptake,'vac_rna2'] <- pop_array[,col_uptake,'vac_rna2'] + uptake_vector   # correct for the negative values
        pop_array[,col_uptake,'vac_adeno2'] <- pop_array[,col_uptake,'vac_adeno2'] - uptake_vector # - negative values == 0
      }
      
      if(any(pop_array[,col_uptake,'vac_rna2']<0)){
        uptake_vector <- pop_array[,col_uptake,'vac_rna2']
        uptake_vector[uptake_vector>0] <- 0 # set all non-negative values to 0
        pop_array[,col_uptake,'vac_booster'] <- pop_array[,col_uptake,'vac_booster'] + uptake_vector
        pop_array[,col_uptake,'vac_rna2'] <- pop_array[,col_uptake,'vac_rna2'] - uptake_vector # - negative values == 0
      }
      
      ### Re-infections ----
      # waning into reinfection path - computed after vac waning to avoid negativity
      new_waning_reinf_adeno1  <- (pop_array[,c('R','Rvoc'),'vac_adeno1'] * 0)
      new_waning_reinf_rna1    <- (pop_array[,c('R','Rvoc'),'vac_rna1'] * 0)  
      new_waning_reinf_adeno2  <- get_new_transitions(pop_array[,c('R','Rvoc'),'vac_adeno2'],ve_waning_infection_booster_rate)  
      new_waning_reinf_rna2    <- get_new_transitions(pop_array[,c('R','Rvoc'),'vac_rna2'],ve_waning_infection_booster_rate)  
      new_waning_reinf_waning  <- get_new_transitions(pop_array[,c('R','Rvoc'),'vac_waning'],ve_waning_infection_booster_rate)  
      new_waning_reinf_booster <- get_new_transitions(pop_array[,c('R','Rvoc'),'vac_booster'],ve_waning_infection_booster_rate)  
      new_waning_reinf_booster_waning    <-get_new_transitions(pop_array[,c('R','Rvoc'),'vac_booster_waning'],ve_waning_infection_booster_rate)  
      
      # loop inside reinfection path
      new_waning_reinf_suppreinfvac    <- get_new_transitions(pop_array[,c('R','Rvoc'),'vac_reinfvac'],ve_waning_infection_rate)  
      

      ##update reinfection
      pop_array[,c('R','Rvoc'),'vac_adeno1'] <- pop_array[,c('R','Rvoc'),'vac_adeno1'] - new_waning_reinf_adeno1 
      pop_array[,'S','vac_reinfvac'] <- pop_array[,'S','vac_reinfvac'] + new_waning_reinf_adeno1[1:10] + new_waning_reinf_adeno1[11:20]
      pop_array[,c('R','Rvoc'),'vac_rna1'] <- pop_array[,c('R','Rvoc'),'vac_rna1']  - new_waning_reinf_rna1
      pop_array[,'S','vac_reinfvac'] <- pop_array[,'S','vac_reinfvac'] + rowSums(new_waning_reinf_rna1)
      pop_array[,c('R','Rvoc'),'vac_adeno2'] <- pop_array[,c('R','Rvoc'),'vac_adeno2'] - new_waning_reinf_adeno2
      pop_array[,'S','vac_reinfvac'] <- pop_array[,'S','vac_reinfvac'] + new_waning_reinf_adeno2[1:10] + new_waning_reinf_adeno2[11:20]
      pop_array[,c('R','Rvoc'),'vac_rna2'] <- pop_array[,c('R','Rvoc'),'vac_rna2'] - new_waning_reinf_rna2
      pop_array[,'S','vac_reinfvac'] <- pop_array[,'S','vac_reinfvac'] + new_waning_reinf_rna2[1:10] + new_waning_reinf_rna2[11:20]
      pop_array[,c('R','Rvoc'),'vac_waning'] <- pop_array[,c('R','Rvoc'),'vac_waning'] - new_waning_reinf_waning
      pop_array[,'S','vac_reinfvac'] <- pop_array[,'S','vac_reinfvac'] + new_waning_reinf_waning[1:10] + new_waning_reinf_waning[11:20]
      pop_array[,c('R','Rvoc'),'vac_booster'] <- pop_array[,c('R','Rvoc'),'vac_booster'] - new_waning_reinf_booster
      pop_array[,'S','vac_reinfvac'] <- pop_array[,'S','vac_reinfvac'] + new_waning_reinf_booster[1:10]  + new_waning_reinf_booster[11:20]
      pop_array[,c('R','Rvoc'),'vac_booster_waning'] <- pop_array[,c('R','Rvoc'),'vac_booster_waning'] - new_waning_reinf_booster_waning
      pop_array[,'S','vac_reinfvac'] <- pop_array[,'S','vac_reinfvac'] + new_waning_reinf_booster_waning[1:10] + new_waning_reinf_booster_waning[11:20] 
      
      pop_array[,col_recovered,'vac_reinfvac'] <- pop_array[,col_recovered,'vac_reinfvac'] - new_waning_reinf_suppreinfvac
      pop_array[,'S','vac_reinfvac'] <- pop_array[,'S','vac_reinfvac'] + new_waning_reinf_suppreinfvac[1:10] + new_waning_reinf_suppreinfvac[11:20] 
      
      }
    
    # Waning immunity after infection ----
    new_waning_reinf                 <- get_new_transitions(pop_array[,c('R','Rvoc'),'orig'],ve_waning_infection_rate)  
    pop_array[,c('R','Rvoc'),'orig'] <- pop_array[,c('R','Rvoc'),'orig'] - matrix(new_waning_reinf,nrow=10)
    pop_array[,'S','vac_reinf']      <- pop_array[,'S','vac_reinf'] + new_waning_reinf[1:10] + new_waning_reinf[11:20]
    new_waning_reinf_suppreinf             <- get_new_transitions(pop_array[,c('R','Rvoc'),'vac_reinf'],ve_waning_infection_booster_rate)
    pop_array[,c('R','Rvoc'),'vac_reinf']  <- pop_array[,c('R','Rvoc'),'vac_reinf'] - matrix(new_waning_reinf_suppreinf,nrow=num_age_groups)
    pop_array[,'S','vac_reinf']            <- pop_array[,'S','vac_reinf'] + new_waning_reinf_suppreinf[1:10] + new_waning_reinf_suppreinf[11:20]
    
    # account for stochastic issues to end up with negative results
    # and if "uptake(a) > alive unvaccinated individuals(a)"
    if(bool_stochastic | bool_vaccine){
      pop_array[pop_array<0] <- 0
    }

    ## Update lists  ---- 
    ##----------------- -
    nfull_array[times_day[step_index+1]+1,,,!pop_categories$vaccine_related] = nfull_array[times_day[step_index+1]+1,,,!pop_categories$vaccine_related] + pop_array[,,!pop_categories$vaccine_related]
    
    # Same for vaccine compartments
    if(bool_vaccine){
      pop_array[pop_array<0] <- 0
      nfull_array[times_day[step_index+1]+1,,,pop_categories$vaccine_related] = nfull_array[times_day[step_index+1]+1,,,pop_categories$vaccine_related] + pop_array[,,pop_categories$vaccine_related]
    }
  }

  ## Model output  ---- 
  ##----------------- - 
  nfull_array_aggr <- get_3dim_nfull_array(nfull_array)
  nsim_day_id <- 0:(ndays_sim-1)
  
  ## Data (log)likelihood  ---- 
  ##------------------------- -

  # aggregate by burden of disease type
  total_new_hosp_icu       <- get_cases_summary(nfull_array_aggr,state_name = 'new_hosp_total')
  total_new_deaths         <- get_cases_summary(nfull_array_aggr,state_name = 'new_D')
  total_nI_sev_sims        <- get_cases_summary(nfull_array_aggr,state_name = 'I_sev')
  total_nHosp_ICU_sims     <- get_cases_summary(nfull_array_aggr,state_name = c('I_hosp','I_icu'))
  total_nHosp_ICU_voc_sims <- get_cases_summary(nfull_array_aggr,state_name = c('Ivoc_hosp','Ivoc_icu'),bool_extend=FALSE)
  total_nICU_sims          <- get_cases_summary(nfull_array_aggr,state_name = 'I_icu')
  total_nHosp_sims         <- get_cases_summary(nfull_array_aggr,state_name = 'I_hosp')
  total_hosp_recoveries    <- get_cases_summary(nfull_array_aggr,state_name = 'hosp_recovered')
  
  #TODO: make collection day flexible
  prob_sero1 <- (serology_sens_function(seq(30,1))%*%nfull_array[1:30,,'asym_mild_infected','orig'])/cohort.size
  prob_sero2 <- (serology_sens_function(seq(51,1))%*%nfull_array[1:51,,'asym_mild_infected','orig'])/cohort.size
  
  total_n1_mild_sims <- rowSums(nfull_array_aggr[,,'Ivoc_mild'])
  total_n2_mild_sims <- rowSums(nfull_array_aggr[,,'I_mild'] + nfull_array_aggr[,,'Ivoc_mild'])
  p_VOC = total_n1_mild_sims/total_n2_mild_sims;
  
  # hospital exit
  total_hosp_exit <- total_hosp_recoveries + total_new_deaths
  
  # fraction in ICU
  p_ICU <- rowSums(total_nICU_sims[,-1]) / rowSums(total_nHosp_ICU_sims[,-1])

  # new infections (i.e. incidence)
  total_new_infections       <- get_cases_summary(nfull_array_aggr,state_name = 'new_E')
  
  # additional statistics
  if(!bool_full_output){
    
    ## Return
    return(list(total_nI_sev_sims = total_nI_sev_sims,
                total_new_infections = total_new_infections,
                total_new_hosp_icu = total_new_hosp_icu,
                total_new_deaths   = total_new_deaths,
                total_nHosp_ICU_sims = total_nHosp_ICU_sims,
                total_nHosp_ICU_voc_sims = total_nHosp_ICU_voc_sims,
                total_nHosp_sims = total_nHosp_sims,
                total_nICU_sims = total_nICU_sims,
                prob_sero1 = prob_sero1,
                prob_sero2 = prob_sero2,
                p_VOC     = p_VOC,
                total_hosp_recoveries = total_hosp_recoveries,
                total_hosp_exit = total_hosp_exit,
                p_ICU=p_ICU))
                
  } else {

    # hospital load
    hosp_icu_load = data.frame(day=nsim_day_id,cases=total_nHosp_ICU_sims[,-1] * h) # mean prevalence per hour
    icu_load      = data.frame(day=nsim_day_id,cases=total_nICU_sims[,-1] * h) # mean prevalence per hour

    # other incidence outcomes
    total_new_reinfections     <- get_cases_summary(nfull_array[,,,pop_categories$reinfection],state_name = 'new_E')
    total_new_mild_infections  <- get_cases_summary(nfull_array_aggr,state_name = 'mild_infected')
    total_recov_mild_infection <- get_cases_summary(nfull_array_aggr,state_name = 'recov_mild_infection')
    total_new_icu              <- get_cases_summary(nfull_array_aggr,state_name = 'new_icu')

    # account for step size per day when inspecting the prevalence
    nfull_array_prev <- nfull_array_aggr[,,] *h
    prev_E           <- get_cases_summary(nfull_array_prev,state_name = 'E')
    prev_I_presym    <- get_cases_summary(nfull_array_prev,state_name = 'I_presym')
    prev_I_asym      <- get_cases_summary(nfull_array_prev,state_name = 'I_asym')
    prev_I_mild      <- get_cases_summary(nfull_array_prev,state_name = 'I_mild')
    prev_I_sev       <- get_cases_summary(nfull_array_prev,state_name = 'I_sev')
    prev_R       <- get_cases_summary(nfull_array_prev,state_name = 'R')
    
    
    # hospital admissions for covid (with estimation/correction)
    pred_hosp1_for_covid      <- as.matrix(total_new_hosp_icu)
    pop_data_be               <- get_regional_pop("belgium")
    pop60                     <- pop_data_be[6] + pop_data_be[7] + pop_data_be[8] + pop_data_be[9] + pop_data_be[10]
    esti_prev60               <- as.matrix(rowSums(prev_I_asym[,7:11])) 
    esti_prev60               <- matrix(esti_prev60, ncol=num_age_groups, nrow=length(esti_prev60))
    pred_hosp1_for_covid[,-1] <- pred_hosp1_for_covid[,-1] - (pred_hosp1_for_covid[,-1]/rowSums(pred_hosp1_for_covid[,-1]))*parms['for_covid_coef']*esti_prev60/pop60
    pred_hosp1_for_covid[pred_hosp1_for_covid < 0] <- 0
    total_new_hosp_icu_for_covid <- pred_hosp1_for_covid
    
    # non-vaccinated hospital admissions
    total_new_hosp_icu_nvac     <- get_cases_summary(nfull_array[,,,!pop_categories$vaccine_related],state_name = 'new_hosp_total')
    
    # vaccinated 2doses and booster hospital admissions
    total_new_hosp_icu_vac_d2     <- get_cases_summary(nfull_array[,,,pop_categories$min_2doses],state_name = 'new_hosp_total')
    
    # vaccinated with booster hospital admissions
    #TODO temporary changed to reinfection admissions for checking !!!!!!!
    total_new_hosp_icu_vac_booster <- get_cases_summary(nfull_array[,,,c('vac_reinf','vac_reinfvac')],state_name = 'new_hosp_total')
    
    # susceptibility states
    mean_susceptible_full               <- get_cases_summary(nfull_array[,,,'orig']*h,state_name = 'S')
    mean_susceptible_vac_d1             <- get_cases_summary(nfull_array[,,,c('vac_rna1','vac_adeno1')]*h,state_name = 'S')
    mean_susceptible_vac_d2             <- get_cases_summary(nfull_array[,,,c('vac_rna2','vac_adeno2')]*h,state_name = 'S')
    mean_susceptible_vac_d2_waning      <- get_cases_summary(nfull_array[,,,'vac_waning']*h,state_name = 'S')
    mean_susceptible_vac_booster        <- get_cases_summary(nfull_array[,,,'vac_booster']*h,state_name = 'S')
    mean_susceptible_vac_booster_waning <- get_cases_summary(nfull_array[,,,'vac_booster_waning']*h,state_name = 'S')
    mean_susceptible_reinf              <- get_cases_summary(nfull_array[,,,'vac_reinf']*h,state_name = 'S')
    mean_susceptible_vac_reinfvac       <- get_cases_summary(nfull_array[,,,'vac_reinfvac']*h,state_name = 'S')
    
    # uptake by day and age
    mean_vaccinated_age  <- get_cases_summary(nfull_array[,,,pop_categories$vaccine_related]*h,state_name = col_no_voc)

    # vaccinated: age, dose type
    i_name <- pop_categories$name[pop_categories$vaccine_related][2]
    for(i_name in pop_categories$name[pop_categories$vaccine_related]){
      vac_age <- get_cases_summary(nfull_array[,,,i_name]*h,state_name = col_no_voc,col_prefix = i_name)
      if(i_name == 'vac_rna1'){
        mean_dose_vaccinated <- vac_age
      } else {
        mean_dose_vaccinated <- cbind(mean_dose_vaccinated,vac_age[,-1])
      }
    }
    names(mean_dose_vaccinated) <- gsub('\\.','_age',names(mean_dose_vaccinated))

    ## Return
    return(list(total_nI_sev_sims = total_nI_sev_sims,
                total_new_hosp_icu = total_new_hosp_icu,
                total_new_icu      = total_new_icu,
                total_new_deaths   = total_new_deaths,
                total_nHosp_ICU_sims = total_nHosp_ICU_sims,
                total_nHosp_ICU_voc_sims = total_nHosp_ICU_voc_sims,
                total_nHosp_sims   = total_nHosp_sims,
                total_nICU_sims = total_nICU_sims,
                prob_sero1 = prob_sero1,
                prob_sero2 = prob_sero2,
                p_VOC      = p_VOC,
                total_hosp_recoveries = total_hosp_recoveries,
                total_hosp_exit = total_hosp_exit,
                p_ICU = p_ICU,
                
                total_new_infections = total_new_infections,
                total_new_mild_infections = total_new_mild_infections,
                total_recov_mild_infection = total_recov_mild_infection,
                hosp_icu_load = hosp_icu_load,
                icu_load = icu_load,
                
                total_new_reinfections = total_new_reinfections,
                
                prev_E=prev_E,
                prev_I_presym=prev_I_presym,
                prev_I_asym=prev_I_asym,
                prev_I_mild=prev_I_mild,
                prev_I_sev=prev_I_sev,
                prev_R=prev_R,
                
                total_new_hosp_icu_nvac = total_new_hosp_icu_nvac,
                total_new_hosp_icu_vac_d2 = total_new_hosp_icu_vac_d2,
                total_new_hosp_icu_vac_booster = total_new_hosp_icu_vac_booster,
                total_new_hosp_icu_for_covid = total_new_hosp_icu_for_covid,
                
                mean_susceptible_full         = mean_susceptible_full,
                mean_susceptible_vac_d1       = mean_susceptible_vac_d1,
                mean_susceptible_vac_d2       = mean_susceptible_vac_d2,
                mean_susceptible_vac_d2_waning= mean_susceptible_vac_d2_waning,
                mean_susceptible_vac_booster  = mean_susceptible_vac_booster,
                mean_susceptible_vac_booster_waning = mean_susceptible_vac_booster_waning,
                mean_susceptible_reinf        = mean_susceptible_reinf,
                mean_susceptible_vac_reinfvac = mean_susceptible_vac_reinfvac,
                
                mean_vaccinated_age = mean_vaccinated_age,
                mean_dose_vaccinated = mean_dose_vaccinated)
    )
  } 
}  

log_likelihood_model <- function(parms, contact_data, method = "mean", plots = "FALSE",
                              parms_names, ndays_sim = NA, #vaccine_uptake = NA,
                              vaccine_uptake = NA,
                              V_mat, V_mat_dose2, V_mat_booster, V_mat_2ndbooster = V_mat_booster*0,  # temp 2nd booster should be added later
                              be_ref_data){

  # run the dynamic model ----
  if(is.na(ndays_sim)) { ndays_sim = nrow(be_ref_data)};
  model_out <- run_model(parms=parms, 
                   contact_data=contact_data, 
                   method=method, 
                   parms_names=parms_names,
                   ndays_sim=ndays_sim,
                   bool_full_output=FALSE,
                   vaccine_uptake = vaccine_uptake,
                   be_ref_data = be_ref_data,
                   V_mat = V_mat,
                   V_mat_dose2 = V_mat_dose2,
                   V_mat_booster = V_mat_booster,
                   V_mat_2ndbooster = V_mat_2ndbooster)
  
    # set number of age groups
    num_age_groups <- 10
  
    # select reference data according to simulation start date
    be_ref_data <- be_ref_data[be_ref_data$date >= get_scm_start_date(),]
    if(min(be_ref_data$date) > get_scm_start_date()){
      date_dummy <- data.frame(date = seq(get_scm_start_date(),max(as.Date(be_ref_data$date)),1))
      be_ref_data <- merge(date_dummy,be_ref_data,by='date',all.x=TRUE)
    }
    
    # get total hospital admissions
    obs_total_hosp_ext   <- be_ref_data$hospital_admissions + be_ref_data$hospital_admissions_other
    obs_total_hosp_ext[is.na(obs_total_hosp_ext)] <- 0
    
    obs_forcovid_hosp    <- be_ref_data$hospital_admissions 
    
    # fix to align previous code
    obs_total_hosp <- obs_total_hosp_ext[-(1:7)]
    
    # age_dist_hosp_mat_ext <- get_regional_hospital_age_distr(be_ref_data$region[1])
    #age_dist_hosp_mat_ext <- get_hospital_age_distr(bool_use_default = FALSE)
    age_dist_hosp_mat_ext <- t(be_ref_data[,grepl('prop_hospital_admission_age',names(be_ref_data))])
    missing_hosp_age_data <- length(obs_total_hosp_ext) - ncol(age_dist_hosp_mat_ext)
    if(missing_hosp_age_data>0){
      age_dist_hosp_mat_ext <- cbind(age_dist_hosp_mat_ext,
                                     age_dist_hosp_mat_ext[,rep(ncol(age_dist_hosp_mat_ext),missing_hosp_age_data)])   
    }
    obs_hosp_load_ext <- be_ref_data$hospital_load
    obs_icu_load_ext  <- be_ref_data$icu_load
    obs_hosp_load_ext[is.na(obs_hosp_load_ext)] <- 0
    obs_icu_load_ext[is.na(obs_icu_load_ext)]   <- 0
    
    # hospital exit
    obs_hosp_exit_ext <- be_ref_data$hospital_exit
    obs_hosp_exit_ext[is.na(obs_hosp_exit_ext)] <- 0
    
    obs_age_dist_mort <- be_ref_data[,grepl('covid19_deaths_age',names(be_ref_data))]
    
    # hospital discharges
    obs_hosp_discharges_ext <- be_ref_data$hospital_discharges
    obs_hosp_discharges_ext[is.na(obs_hosp_discharges_ext)] <- 0
    
    # hospital recoveries: estimated by load, admissions and mortality
    obs_hosp_recoveries_ext <- be_ref_data$hospital_recovery_7dmean
    obs_hosp_recoveries_ext[is.na(obs_hosp_recoveries_ext)] <- 0
    
    # Serological survey data
    sero_data_all <- read.table('data/covid19_data_serology_20230526.csv',sep=',',header=T)
     
    # VOC data
    db_voc_raw <- read.table('data/covid19_data_voc_20230526.csv',sep=',',header=T,stringsAsFactors = T)
    db_voc_raw$date <- as.Date(db_voc_raw$date)
    
    # end: select reference data
  
  # re-capture some parameters 
  h        = parms[parms_names == 'h'];  
  sel_region = get_region(parms[parms_names == 'region_id'])
  
  # parse model output
  total_nI_sev_sims      = model_out[['total_nI_sev_sims']]
  total_new_hosp_icu     = model_out[['total_new_hosp_icu']]
  total_new_hosp_icu_for_covid     = model_out[['total_new_hosp_icu_for_covid']]
  total_new_deaths       = model_out[['total_new_deaths']]
  total_nHosp_ICU_sims   = model_out[['total_nHosp_ICU_sims']]       # required for Binomial distribution
  total_nHosp_ICU_voc_sims = model_out[['total_nHosp_ICU_voc_sims']] # required for Binomial distribution
  total_nHosp_sims       = model_out[['total_nHosp_sims']]           # required for Binomial distribution
  total_nICU_sims        = model_out[['total_nICU_sims']]            # required for Binomial distribution
  
  total_nHosp_ICU_load = rowSums(total_nHosp_ICU_sims[,-1])*h  # = average hospital load / day
  total_nICU_load    = rowSums(total_nICU_sims[,-1])*h  # = average ICU load / day
  prob_sero1         = model_out[['prob_sero1']]
  prob_sero2         = model_out[['prob_sero2']]
  p_VOC             = model_out[['p_VOC']]

  total_hosp_recoveries = rowSums(model_out[['total_hosp_recoveries']][,-1])
  total_hosp_exit       = rowSums(model_out[['total_hosp_exit']][,-1])
  p_ICU                 = model_out[['p_ICU']]
  
  omega    = exp(parms[grepl('log_omega',parms_names)]); 
  delta3   = exp(parms[grepl('log_delta3',parms_names)]); 
  delta3_wave2 = exp(parms['log_delta3_wave2']);    # recovery rate from hospital (from sept 2021)
  
  total_new_infections <- model_out[['total_new_infections']]
  cumsum_new_infections_age <- apply(total_new_infections[,-1],2,cumsum)
  
  ### 1. Data sources (hospitalizations, mortality, serology)  ---- 
  ###------------------------------------------------------------ -
  
  cohort.size <- get_regional_pop(sel_region)
  
  # ## REFERENCE DATA
  # be_ref_data <- get_latest_incidence_data()
  
  #### 1.1. Daily number of new hospitalizations by age group (from 1/3 onwards)
  ####------------------------------------------------------------------------ -
  y1 = obs_total_hosp_ext;
  y1_mat = round_tot(y1*t(age_dist_hosp_mat_ext),0)
  if(nrow(y1_mat)>ndays_sim){
    y1_mat <- y1_mat[1:ndays_sim,]
  }
  
  #### 1.2. Daily number of new deaths by age group (from 1/3 onwards) 
  ####-------------------------------------------------------------- -
  y2_mat = obs_age_dist_mort;
  y2 = rowSums(y2_mat)
  if(nrow(y2_mat)>ndays_sim){
    y2_mat <- y2_mat[1:ndays_sim,]
  }
  
  #### 1.3. Serial serological survey results based on residual samples 
  ####--------------------------------------------------------------- -
  y3 = sero_data_all$igg_cat_pos[sero_data_all$cround == 1];
  n3 = sero_data_all$igg_cat_total[sero_data_all$cround == 1];
  y4 = sero_data_all$igg_cat_pos[sero_data_all$cround == 2];
  n4 = sero_data_all$igg_cat_total[sero_data_all$cround == 2];
  
  # using the observed/modelled incidence data
  y3_infection      <- apply(be_ref_data[,grepl('cases_',names(be_ref_data))],2,cumsum)
  n3_infection      <- cohort.size
  y3_infection_days <- 1:min(365,                 # one year, to prevent issues with re-infections
                             nrow(y3_infection),  # available reference data
                             ndays_sim)           # simulation time horizon # until first missing value
  
  # restrict time horizon
  y3_infection              <- y3_infection[y3_infection_days,]
  cumsum_new_infections_age <- cumsum_new_infections_age[y3_infection_days,]
  
  #### 1.4. Prevalence of Variants of Concern: VOC strain in the model (alpha, omicron)
  ####----------------------------------------- -
  y5 = db_voc_raw$voc_aggr + db_voc_raw$voc_omicron
  n5 = db_voc_raw$n_sequenced
  y5_date <- db_voc_raw$date
  y5_days <- sim_date2day(y5_date)
  
  # make daily
  y5_days <- sim_date2day(approx_sequenced(y5_date,y5)$x)
  y5      <- approx_sequenced(y5_date,y5)$y
  n5      <- approx_sequenced(y5_date,n5)$y

  # adjust nday_sim
  if(max(y5_days)>ndays_sim){
    y5     <- y5[y5_days<=ndays_sim]
    n5     <- n5[y5_days<=ndays_sim]
    y5_days <- y5_days[y5_days<=ndays_sim] 
  }
  
  ####----------------------------------------- -
  y6_exit <- obs_hosp_exit_ext
  if(length(y6_exit)>ndays_sim){
    y6_exit <- y6_exit[1:ndays_sim]
  }
  
  y6_load <- obs_hosp_load_ext
  if(length(y6_load)>ndays_sim){
    y6_load <- y6_load[1:ndays_sim]
  }
  
  #### 1.6. Hospital and ICU load
  ####----------------------------------------- -
  y7_icu  <- obs_icu_load_ext  
  n7_hosp <- obs_hosp_load_ext 
  y7_days <- which(!is.na(y7_icu))
  if(length(y7_icu)>ndays_sim){
    y7_icu  <- y7_icu[1:ndays_sim]
    n7_hosp <- n7_hosp[1:ndays_sim]
    y7_days <- y7_days[y7_days<=ndays_sim] 
  }
  
  #### 1.1. Daily number of new hospitalizations for covid
  ####------------------------------------------------------------------------ -
  obs_forcovid_hosp <- obs_forcovid_hosp[1:ndays_sim]
  
  ### 2. Loglikelihood contributions  ---- 
  ###----------------------------------- -
  
  #### 2.1. Binomial distribution for the daily number of new hospitalizations 
  ####---------------------------------------------------------------------- -
  ll_hosp_icu_mat <- matrix(nrow = nrow(y1_mat), ncol = num_age_groups)
  for (j in 1:num_age_groups){
      ll_hosp_icu_mat[,j] = dbinom(y1_mat[,j], size = round(unlist(total_nI_sev_sims[1:nrow(y1_mat),j+1]),0), 
                                 prob = 1 - exp(-h*omega[j]), log = T)  
  }
  
  ll_hosp_icu_mat = ifelse(ll_hosp_icu_mat == -Inf, -100000, ll_hosp_icu_mat)
  
  #### 2.2. Weighed Least Squares for the daily number of new deaths (assuming equal death rates for hosp/ICU)
  ####------------------------------------------------------------------------------------------------------ -
  #note: transition rates change over time with adjusted LoS, mortality, severity...
 
  ls_hosp_mort  = sum((y2_mat - total_new_deaths[1:nrow(y2_mat),-1])^2,na.rm=T)
  wls_hosp_mort = sum((y2_mat[,-1] - total_new_deaths[1:nrow(y2_mat),-(1:2)])^2 / total_new_deaths[1:nrow(y2_mat),-(1:2)],na.rm=T)
  
  ll_mort_mat    <- matrix(nrow = nrow(y2_mat), ncol = num_age_groups)
  j <- 2
  for (j in 2:num_age_groups){
    if(sum(y2_mat[,j],na.rm=T)/(ndays_sim/7)<1){
      rollmean_k <- 30
    } else{
      rollmean_k <- 7
    }
    n_mort    <- round(rollsum(unlist(total_new_deaths[1:nrow(y2_mat),j+1]),k=rollmean_k,na.rm=T,fill = NA),0)
    hosp_load <- round(rollsum(unlist(total_nHosp_ICU_sims[1:nrow(y2_mat),j+1]),k=rollmean_k,na.rm=T,fill = NA),0)
    p_mort    <- n_mort / hosp_load
    
    ll_mort_mat[,j] <- dbinom(x = rollsum(y2_mat[,j],k=rollmean_k,na.rm=T,fill = NA),
                              size = hosp_load, 
                              prob = p_mort, 
                              log = T)
  }
  
  ll_mort_mat = ifelse(ll_mort_mat == -Inf, -100, ll_mort_mat)

  #### 2.3.a Binomial distribution for the serial serological data 
  ####--------------------------------------------------------- -
  ll_sero1 = dbinom(y3, size = n3, prob = prob_sero1, log = T)
  ll_sero1 = ifelse(ll_sero1 == -Inf, -100000, ll_sero1)
  ll_sero2 = dbinom(y4, size = n4, prob = prob_sero2, log = T)
  ll_sero2 = ifelse(ll_sero2 == -Inf, -100000, ll_sero2)
  
  #### 2.3.b Binomial distribution based on (reported/simulated) new infections 
  ####--------------------------------------------------------- -
  ll_infections <- matrix(nrow = length(y3_infection_days), ncol = num_age_groups)
  for(i_age in 1:num_age_groups){
    ll_infections[,i_age] = dbinom(y3_infection[,i_age], size = n3_infection[i_age], prob = cumsum_new_infections_age[,i_age]/cohort.size[i_age], log = T)
  }
  ll_infections = ifelse(ll_infections == -Inf, -100000, ll_infections)
  
  #### 2.4. Binomial distribution for the VOC-strain prevalence (start 1/1/2021)
  ####------------------------------------------------------------------ -
  ll_VOC = dbinom(y5, size = n5, prob = p_VOC[y5_days], log = T);
  ll_VOC = ll_VOC[ll_VOC != -Inf]
  
  
  #### 2.4.a Binomial distribution for the hospital exit
  ####------------------------------------------------------------------ -
  ll_hosp_exit = dbinom(y6_exit, size = round(rowSums(total_nHosp_ICU_sims[1:length(y6_exit),-1]),0), 
                             prob = 1 - exp(-h*delta3), log = T)
  ll_hosp_exit = ifelse(ll_hosp_exit == -Inf, -100000, ll_hosp_exit)
  
  #### 2.4.b Weighed Least Squares for the hospital load
  ####------------------------------------------------------------------ -
  ls_hosp_load  = sum((y6_load - total_nHosp_ICU_load[1:length(y6_load)])^2,na.rm=T)
  wls_hosp_load = sum((y6_load - total_nHosp_ICU_load[1:length(y6_load)])^2 / total_nHosp_ICU_load[1:length(y6_load)],na.rm=T)
  
  #### 2.5. Binomial distribution for the fraction ICU load / hospital load
  ####------------------------------------------------------------------ -
  ll_icu_fraction = dbinom(y7_icu[y7_days], size = n7_hosp[y7_days], prob = p_ICU[y7_days], log = T)
  ll_icu_fraction = ifelse(ll_icu_fraction == -Inf, -100000, ll_icu_fraction)
  
  #### 2.4.b Weighed Least Squares for the ICU load
  ####------------------------------------------------------------------ -
  wls_icu_load = sum((y7_icu[y7_days] - total_nICU_load[y7_days])^2 / total_nICU_load[y7_days],na.rm=T)
  
  
  #### Weighed Least Squares for "for covid" coef
  ####------------------------------------------------------------------ -
  ls_forcovid  = sum((obs_forcovid_hosp - total_new_hosp_icu_for_covid[1:length(obs_forcovid_hosp)])^2,na.rm=T)
  wls_forcovid = sum((obs_forcovid_hosp - total_new_hosp_icu_for_covid[1:length(obs_forcovid_hosp)])^2 / total_new_hosp_icu_for_covid[1:length(obs_forcovid_hosp)],na.rm=T)
  
  
  ## Graphical exploration  ---- 
  ##-------------------------- - 
  if(plots == "TRUE"){
    plot_log_likelihood_model(total_nI_sev_sims,total_new_hosp_icu,
                           total_new_deaths,total_nHosp_ICU_load,
                           total_nICU_load,
                           total_hosp_exit,
                           y1_mat,y2_mat,
                           prob_sero1,prob_sero2,y3,n3,y4,n4,
                           y5,n5,y5_days,p_VOC,n7_hosp,y7_icu,
                           y6_exit,y6_days,
                           cumsum_new_infections_age,
                           y3_infection,n3_infection)
  }
  
  ## Loglikelihood function value  ---- 
  ##--------------------------------- -
  w1 = (cohort.size/n3);         # poststratification weights (inverse probability weighing)
  w2 = (cohort.size/n4);         # poststratification weights (inverse probability weighing)
  
  w1_imp = 1                                                      # importance weights
  w3_imp = 1/(length(ll_sero1)/length(ll_hosp_icu_mat))           # importance weights
  w4_imp = 1/(length(ll_sero2)/length(ll_hosp_icu_mat))           # importance weights
 
  ### Weighted likelihood function: inverse probability and importance weighting
  ###----------------------------------------------------------- -
  sumll1 = w1_imp*sum(ll_hosp_icu_mat,na.rm=T) + 
           w3_imp*sum((w1/sum(w1))*ll_sero1,na.rm=T) + 
           w4_imp*sum((w2/sum(w2))*ll_sero2,na.rm=T);
  
  ### Weighted likelihood function: hosp+sero+VOC
  ###-------------------------------------------------------------------------- -
  sumll2 = sumll1 + w1_imp*sum(ll_VOC,na.rm=T); 
  
  ### Weighted likelihood function: hospital and ICU load
  ###-------------------------------------------------------------------------- -
  sumll3 = - wls_hosp_load - wls_icu_load
  
  ### WLS function: mortality
  ###-------------------------------------------------------------------------- -
  sumll4 = sum(ll_mort_mat,na.rm=T)
  
  
  ### Weighted likelihood function: VOC
  ###-------------------------------------------------------------------------- -
  sumll5 = sum(ll_VOC,na.rm=T)
  
  ### Weighted likelihood function: for covid + load + mortality
  ###-------------------------------------------------------------------------- -
  sumll6 = -wls_forcovid + sumll3 + sumll4
  
  ### likelihood function: infections and hospitalisation
  ###-------------------------------------------------------------------------- -
  sumll7 = sum(ll_hosp_icu_mat,na.rm=T) + 
            sum(ll_infections,na.rm=T)
  
  return(list(crit1 = sumll1, 
              crit2 = sumll2, 
              crit3 = sumll3, 
              crit4 = sumll4, 
              crit5 = sumll5,
              crit6 = sumll6,
              crit7 = sumll7,
              dev1 = -2*sumll1,
              dev2 = -2*sumll2,
              dev3 = -2*sumll3,
              dev4 = -2*sumll4,
              dev5 = -2*sumll5,
              dev6 = -2*sumll6,
              dev7 = -2*sumll7
              ))
}

## Graphical exploration  ---- 
##-------------------------- - 
plot_log_likelihood_model <- function(total_nI_sev_sims,total_new_hosp_icu,
                                   total_new_deaths,total_nHosp_ICU_load,
                                   total_nICU_load,
                                   total_hosp_exit,
                                   y1_mat,y2_mat,
                                   prob_sero1,prob_sero2,y3,n3,y4,n4,
                                   y5,n5,y5_days,p_VOC,n7_hosp,y7_icu,
                                   y6_exit,y6_days,
                                   cumsum_new_infections_age,
                                   y3_infection,n3_infection,
                                   bool_seroprev=FALSE){

      # check figure margins, and abort if too large (or figure is to narrow)
      if(any(dev.size()<4)){
        warning('Issue with dev.size, figure are too narrow,, skip plotting')
        return(NULL)
      }
  
      if(!is_vsc()) par(mfrow = c(3,3))
  
      plot_model_output <- function(run_data_age,
                                    ref_data_age,
                                    y_lab = 'New cases',
                                    x_lab = "Days since introduction",
                                    bool_by_age = TRUE){
        
        # make sure we work with data.frames (hence a matrix on which we can call rowSums)
        run_data_age <- data.frame(run_data_age)
        ref_data_age <- data.frame(ref_data_age)
        
        if(ncol(run_data_age)!= ncol(ref_data_age)){
          warning('plot_model_output(): issue with dimensions of run and reference data')
        }
        
        # aggregated
        x_lim <- range(0,nrow(ref_data_age),nrow(run_data_age))
        y_lim <- c(0,max(rowSums(ref_data_age,na.rm = T),
                          rowSums(run_data_age,na.rm = T))*1.1)
        plot(1:nrow(ref_data_age), rowSums(ref_data_age),cex=0.5,
             ylab = y_lab, 
             xlab = x_lab,
             ylim = y_lim,
             xlim = x_lim)
        lines(rowSums(run_data_age), col = 2, lwd = 2)
        legend('topleft',
               c('reference','model'),
               pch=c(1,NA),
               lwd=c(NA,2),
               col = 1:2,
               cex=0.7)
        
        # by age
        if(bool_by_age && ncol(run_data_age) == 10){
          y_lim <- c(0,max(c(unlist(ref_data_age)),
                           c(unlist(run_data_age)),na.rm=T)*2)
          p_age <- seq(1,ncol(ref_data_age),2)
          p_col = rep(1:length(p_age),each=2)
          plot(0, 0,cex=0.5,
               ylab = y_lab, 
               xlab = x_lab,
               col  = 0,
               ylim = y_lim,
               xlim = x_lim)
          for(i_age in p_age){
            col_age <- i_age + (0:1)
            points(1:nrow(ref_data_age), rowSums(ref_data_age[,col_age]),cex=0.5,pch=i_age,col=p_col[i_age])
            lines(rowSums(run_data_age[,col_age]), col = p_col[i_age], lwd = 2)
          }
          legend('topleft',
                 paste(p_age,p_age+1,sep='+'),
                 pch=16,
                 ncol=5,
                 cex=0.7,
                 col = p_col[p_age]) 
        }
      }
  
      # hospital admissions: total and by age
      plot_model_output(run_data_age = total_new_hosp_icu[,-1],
                        ref_data_age = y1_mat,
                        y_lab = "Number of new hospitalizations")
      
      
      
      # mortality: total
      plot_model_output(run_data_age = total_new_deaths[,-1],
                        ref_data_age = y2_mat,
                        y_lab = "Number of new deaths",
                        bool_by_age=FALSE)
      
      # VOC
      p_VOC_ref <- rep(NA,length(p_VOC))
      p_VOC_ref[y5_days] <- y5/n5
      plot_model_output(run_data_age = p_VOC,
                        ref_data_age = p_VOC_ref,
                        y_lab = "Prevalence of VOC-strain")
      
      ## HOSPITAL LOAD
      plot_model_output(run_data_age = total_nHosp_ICU_load,
                        ref_data_age = n7_hosp,
                        y_lab = "Hospital load")
      
      ## ICU LOAD
      plot_model_output(run_data_age = total_nICU_load,
                        ref_data_age = y7_icu,
                        y_lab = "Hospital load")
      
      ## INFECTIONS: total and by age
      plot_model_output(run_data_age = cumsum_new_infections_age,
                        ref_data_age = y3_infection,
                        y_lab = "Cumulative number of infections",
                        x_lab = "Days since introduction\n[YEAR 1]")
      
      ## SEROPREVALENCE
        plot(1:length(prob_sero1),y3/n3, ylim = c(0,0.2),
             xlab = "Age category", ylab = "Seroprevalence",pch=16)
        lines(1:length(prob_sero1), prob_sero1, col = 1, lwd = 2)
        
        points(1:length(prob_sero2),y4/n4,col=4,pch=16)
        lines(1:length(prob_sero2), prob_sero2, col = 4, lwd = 2)
        legend('topleft',
               c('March 2020',
                 'April 2020'),
               pch=1,
               lwd=2,
               cex=0.7,
               col=c(1,4),
               ncol=2)
}


## Log prior distributions
##------------------------ -
log_prior_model <- function(parms,parms_names){
  
  # make sure that the parameter names are set
  names(parms[1:length(parms_names)]) <- parms_names
  
  # set number of age groups
  num_age_groups <- 10
  
  gamma  = exp(parms['log_gamma']);     # average length of the latency period - 2 days
  theta  = exp(parms['log_theta']);     # infectious period at pre-symptomatic stage - 3.5 days
  delta1 = exp(parms['log_delta1']);    # infectious period at asymptomatic stage - 3.5 days
  delta2 = exp(parms['log_delta2']);    # infectious period at (mild) symptomatic stage - 3.5 days
  
  phi0  = expit(parms[grepl('log_phi0',parms_names)]); # proportion of symptomatically infected with mild symptoms
  p_vec = parms[grepl('p_asympt_age',parms_names)];    # proportion of asymptomatic cases
  
  
  omega = exp(parms[grepl('log_omega',parms_names)]);              # waiting time between symptom onset and hospitalization
  phi1  = expit(parms[grepl('log_phi1_age',parms_names)]);         # proportion of severely infected with regular hospitalization

  nStage      <- sum(grepl('log_q_stage',parms_names)) / num_age_groups 
  q_all       <- exp(parms[grepl('log_q_stage',parms_names)]); 
  q_stage     <- matrix(q_all, nrow = nStage, ncol = num_age_groups, byrow = T);

  # all mu parameters
  mu_param = parms[grepl('mu.*_sev',parms_names)]
  mu_sev = expit(mu_param);          # age-specific hospital fatality ratio - probability of dying
  delta3 = exp(parms[grepl('log_delta3',parms_names)]);  # recovery period in severely infected/length of stay - 7 days

  mortality_changepoint = parms[grepl('mortality_changepoint',parms_names)] #FYI: not part of estimation process
  
  n0 = exp(parms["log_n0"]) # imported cases in 2020
  
  voc_start = parms[grepl('VOC.*_start',parms_names)]
  
  voc_init    = exp(parms[grepl('log.*_init',parms_names)]) #exp(parms['log_VOC_alpha_init']);
  voc_transm  = exp(parms[grepl('log.*_transm',parms_names)]) #exp(parms['log_VOC_alpha_transm']);    # additional transmissibility
  voc_hr_hosp = exp(parms[grepl('log.*hr_hosp',parms_names)]) # hospitalization admission hazard ratio
  
  phi1_add  = parms[grepl('phi1_add',parms_names)]  # additional factor for regular hospitalization (instead of ICU)
  
  voc_gamma_factor = exp(parms[grepl('log.*_VOC_.*_gamma_factor',parms_names)]); # adjustment factor for the latency period with Omicron
  
  log_gamma_prior    <- dnorm(gamma, mean = 1/2, sd = 0.05, log = T);
  log_theta_prior    <- dnorm(theta, mean = 1/3.5, sd = 0.05, log = T);
  log_delta1_prior   <- dnorm(delta1, mean = 1/3.5, sd = 0.05, log = T);
  log_delta2_prior   <- dnorm(delta2, mean = 1/7, sd = 0.05, log = T);
  
  log_phi0_prior     <- dunif(phi0, 0.0000001, 0.99999, log = T);
  log_omega_prior    <- dunif(omega, 1/7 , 1/0.5, log= T)
  log_phi1_prior     <- dunif(phi1, 0.0000001, 0.99999, log = T);
  
  log_mu_prior       <- dunif(mu_sev, 0.0000001, 0.99999, log = T);
  log_delta3_prior   <- dunif(delta3, 0, 0.2, log = T);
  log_phi1_add_prior <- dunif(phi1_add, 0, 3, log = T);
  log_n0_prior       <- dunif(n0, 0, 5000, log = T);
  
  log_voc_start        <- dunif(voc_start, global_lib_voc$start_ll[1], global_lib_voc$start_ul[1], log = T);
  for(i_voc in 1:length(voc_start)){
    log_voc_start[i_voc]       <- dunif(voc_start[i_voc], global_lib_voc$start_ll[i_voc], global_lib_voc$start_ul[i_voc], log = T);
  }
  
  log_voc_init_prior   <- dunif(voc_init, 20, 10000, log = T);
  log_voc_transm_prior <- dunif(voc_transm, 0.001, 11, log = T);
  log_voc_hr_hosp_prior<- dunif(voc_hr_hosp, 0, 15, log = T);
  
  log_voc_gamma      <- dunif(voc_gamma_factor, 1,1e10, log = T)
  
  # social contact data
  contact_mat <- readRDS('data/social_contact_data.rds')
  
  fparam    <- parms['f_asympt']
  betas_asy <- fparam * q_stage[1,]  * select_contact_matrix(contact_mat,1,bool_asymptomatic = TRUE)
  betas_sy  <- q_stage[1,]  * select_contact_matrix(contact_mat,1,bool_asymptomatic = FALSE)
  
  p_asy <- p_vec;        
  D1    <- 1/theta; 
  D2    <- 1/delta1;
  D3    <- 1/delta2;
  D4    <- 1/omega;
  
  # R0 estimation is based on national statistics
  cohort.size <- get_regional_pop(region = "belgium")
  
  R1 <- D1*(betas_asy * cohort.size) + 
        p_asy*D2*(betas_asy * cohort.size) + 
        (1-p_asy)*D3*(betas_sy * cohort.size) +
        (1-p_asy)*(1-phi0)*D4*(betas_sy * cohort.size);
  
  R0            <- Re((eigen(R1)$values[1]));
  
  log_q_prior        <- dunif(q_all, exp(-7.7), exp(0.27), log = T);
  
  log_prior_val <-  c(log_gamma_prior, 
                    log_theta_prior, 
                    log_delta1_prior, 
                    log_delta2_prior,
                    log_q_prior, 
                    log_phi0_prior, 
                    log_omega_prior,
                    log_mu_prior,
                    log_delta3_prior, 
                    #log_alpha_prior,
    
                    log_n0_prior,
    
                    log_voc_start,
                    
                    log_voc_init_prior,
                    log_voc_transm_prior, 
                    log_voc_hr_hosp_prior,
                    log_phi1_add_prior,
                    
                    log_voc_gamma)
  
  ll_penalty <- -500000
  log_prior_val <- ifelse(log_prior_val == -Inf, ll_penalty, log_prior_val)

  invalid_param <- NULL
  # include check  
  if(any(log_prior_val == ll_penalty)){
    invalid_param <- names(log_prior_val)[log_prior_val == ll_penalty]
  }

  log_prior_val <- sum(log_prior_val)
  
  return(list(log_prior_val = log_prior_val, R0val = R0, invalid_param = invalid_param))
}


