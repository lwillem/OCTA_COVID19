########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains the main function to run multiple simulation runs with 
# different social contact and/or vaccination uptake assumptions for Belgium 
# or one of the three regions. 
#
# Current parallelism: for different vaccine uptake file
#
# Copyright 2024, SIMID                             
########################################################################### #
#db_cnt_adjust = NULL; adjusted_parameters = NULL;cl_args = NULL; sel_region = NULL;contact_data = NA;be_ref_data = NA;vacc_files = NA
# define function
projections_vaccination <- function(output_tag,
                                    sel_region = NULL,
                                    chains_param_files,
                                    db_cnt_adjust = NULL,
                                    scen_tag,
                                    num_param_sets,
                                    num_stochastic_real,
                                    num_days_sim,
                                    vacc_files = NA,
                                    adjusted_parameters = NULL,
                                    cl_args = NULL,
                                    contact_data = NA,
                                    be_ref_data = NA){

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# use function arguments

# check chains_param_files
if(length(chains_param_files) == 0){
  stop('ERROR: no parameter file provided, stop model.')
}    
# use default y axis?
default_y_axis <- TRUE

################################################################ #
# LOAD AND PRE-PROCESS ----
################################################################ #

# adjust parameters based on optional command line arguments
if(length(cl_args)>=1){ # region
  sel_region <- cl_args[1]
  }
 
# optional: set default region
if(is.null(sel_region)){
  sel_region <- get_region(1)
}

# include contact scenario info into output tag. If NULL, leave empty
cnt_adjust_tag <- paste0(round(db_cnt_adjust$value*100),collapse='_')
if(nchar(cnt_adjust_tag)>0){
  cnt_adjust_tag <-paste0('_cnt',cnt_adjust_tag)
}
if(!all(is.na(contact_data)) && !is.null(contact_data$label)){
  cnt_adjust_tag <- paste0('_',contact_data$label)
}

# select and load parameter file
chains_param_file <- chains_param_files[grepl(sel_region,chains_param_files) ]

if(length(chains_param_file) == 0){
  chains_param_file <- chains_param_files[1] 
}

print(chains_param_file)
parms_chains     = read.table(chains_param_file, sep = ",", header = T)

# select the parameter sets and make a copy for each stochastic realisations
parms_chains  <- parms_chains[rep(1:num_param_sets,num_stochastic_real),]
num_runs      <- nrow(parms_chains)

# option to adjust some specific parameters
if(!is.null(adjusted_parameters) && length(adjusted_parameters)>0){
  for(i_param in 1:length(adjusted_parameters)){
    if(names(adjusted_parameters)[i_param] %in% names (parms_chains)){
      parms_chains[names(adjusted_parameters)[i_param]] <- adjusted_parameters[i_param]
    } else{
      warning('Issue in function projections_vaccination: "adjusted_parameters" contains a parameter that is not present in "parms_chains"')
    }
  }
}

# set parameter names
parms_names <- names(parms_chains)
if(!exists('sel_region')){
  sel_region <- get_region(unique(parms_chains[,'region_id']))
} else{
  parms_chains[,'region_id'] <- get_region_id(sel_region)
}

# set output directory
if(grepl('data/config/',chains_param_file)){
  exp_dir <- gsub('.*p.*crit._','',chains_param_file)
  exp_dir <- gsub('.csv','',exp_dir)
  exp_dir <- basename(exp_dir)
} else{
  exp_dir <- basename((dirname(chains_param_file)))
}

# include experiment details in directory name
exp_dir <- paste0(exp_dir,'/',scen_tag,'_d',num_days_sim,cnt_adjust_tag,'_c',num_param_sets,'_s',num_stochastic_real,'_n',num_runs,'_',sel_region)

# get output directory (and initiate if needed)
output_dir <- init_output_dir(output_tag,exp_dir)
print(output_dir)

# add adjusted contact behaviour to the model parameter matrix
if(!is.null(db_cnt_adjust)){
  for(i_cp in 1:nrow(db_cnt_adjust)){
    if(db_cnt_adjust$value[i_cp]!=0){
      parms_chains[,paste0('cnt_adjust_value',i_cp)]       <- db_cnt_adjust$value[i_cp]
      parms_chains[,paste0('cnt_adjust_day',i_cp)]         <- sim_date2day(db_cnt_adjust$date)[i_cp]
      parms_chains[,paste0('cnt_adjust_reference_day',i_cp)]   <- sim_date2day(db_cnt_adjust$reference_date)[i_cp]
    }
  }
}

# add adjusted parameter config to output folder
write.table(parms_chains,
            file=file.path(output_dir,basename(chains_param_file)),
            sep=',',row.names = F)


# load reference data, if not given as function argument
if(all(is.na(be_ref_data))){
  be_ref_data <- get_latest_incidence_data(sel_region = sel_region)
}

# load social contact data
if(any(is.na(contact_data))){
  contact_data <- readRDS('data/social_contact_data.rds')
}

# make sure the q_stage_id is numeric
contact_data$db_C_sim$q_stage_id <- as.numeric(contact_data$db_C_sim$q_stage_id)

# check compatibility of contact data and q parameters
if(identify_total_nb_stages(names(parms_chains)) < max(contact_data$db_C_sim$q_stage_id)){
  stop('Number of q-parameters does not align with provided contact data')
}

## select and load vaccine uptake
vacc_files <- vacc_files[(grepl(sel_region,vacc_files) )]
vacc_file_tag <- tolower(gsub('\\.csv','',gsub('.*Output_v','',vacc_files)))
vacc_file_tag <- basename(vacc_file_tag) #tmp fix when using different uptake files
if(length(vacc_files)==0){
  vacc_files <- NA
  vacc_file_tag <- 'novac'
}
if(!any(is.na(vacc_files))){
  print(vacc_files)
}

# set default plot color and pch
col_lib <- data.frame(tags = c('xxxx',
                               cnt_adjust_tag,
                               'mild_alpha',
                               'mild_delta'),
                      col  = c('black','darkred','darkblue','darkgreen'),
                      label = c(paste0('Reported (',max(be_ref_data$date),')'),
                                paste('Scenario:',cnt_adjust_tag),
                                'Alpha+Beta+Gamma',
                                'Delta'),
                      lwd = c(NA,2,2,2),
                      pch = c(1,NA,NA,NA))

# set colors for post-processing figures
vacc_colors    <- vacc_file_tag
is_invalid_col <- !vacc_colors %in% colors()
vacc_colors[is_invalid_col] <- (1:sum(is_invalid_col))+1

################################################################ #
# MAIN: CALL MODEL ----
################################################################ #

# explore contact and Q(age) parameters
pdf(file=file.path(output_dir,paste0('q_parameters',gsub('\\.','_p',cnt_adjust_tag),'.pdf')),14,7)
explore_q_param(parms_chains = parms_chains,
                 db_C_sim     = contact_data$db_C_sim)
dev.off()

# # start parallel nodes ----
smd_start_cluster(num_proc = min(length(vacc_files),6),
                  timeout  = 600) # 10min

# run ----
i_vac <- 1
foreach(i_vac = 1:length(vacc_files),
        .export = c('par_nodes_info'),
        .packages = c('EpiEstim',
                      'simid.rtools',
                      'zoo',
                      'openxlsx',
                      'scales',
                      'RColorBrewer')
        ) %dopar%{
  
  # temp: to cope with the parallel environment
  source('R/main_vaccination.R')

  i_vac_file        <- vacc_files[i_vac]
  if(!is.na(i_vac_file)){
    vacc_schedule <- read.csv(i_vac_file, header = T)
  } else{
    vacc_schedule <- NA
  }

  vacc_schedule_tag <- vacc_file_tag[i_vac]
  col_scen          <- vacc_colors[i_vac]
  
  # fix for color_"tag"
  col_scen          <- unlist(strsplit(col_scen,'_'))[1]
  
  # initialise summary variables
  i_seq             <- round(seq(1,nrow(parms_chains),length.out = num_runs))
  hosp_age_adm      <- array(dim = c(num_runs, num_days_sim, 11))
  hosp_age_adforcovid      <- array(dim = c(num_runs, num_days_sim, 11))
  scen_cases        <- array(dim = c(num_runs, num_days_sim, 11))
  scen_mild_cases   <- array(dim = c(num_runs, num_days_sim, 11))
  mortality_age     <- array(dim = c(num_runs, num_days_sim, 11))
  cumul_cases       <- array(dim = c(num_runs, num_days_sim, 11))
  
  prev_E_sims       <- array(dim = c(num_runs, num_days_sim, 10)) 
  prev_I_presymp_sims <- array(dim = c(num_runs, num_days_sim, 10)) 
  prev_I_asymp_sims <- array(dim = c(num_runs, num_days_sim, 10))
  prev_I_mild_sims  <- array(dim = c(num_runs, num_days_sim, 10))
  prev_I_severe_sims<- array(dim = c(num_runs, num_days_sim, 10))
  prev_I_hosp_sims  <- array(dim = c(num_runs, num_days_sim, 10))
  prev_I_icu_sims   <- array(dim = c(num_runs, num_days_sim, 10)) 
  prev_R_sims       <- array(dim = c(num_runs, num_days_sim, 10)) 
  prev_M_sims       <- array(dim = c(num_runs, num_days_sim, 10)) 
  inc_I_hosp_sims   <- array(dim = c(num_runs, num_days_sim, 10))
  inc_I_icu_sims    <- array(dim = c(num_runs, num_days_sim, 10)) 
  inc_I_hosp_for_covid_sims   <- array(dim = c(num_runs, num_days_sim, 10))
  
  mean_susceptible_full         <- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_vac_d1       <- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_vac_d2       <- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_vac_d2_waning<- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_vac_booster  <- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_vac_booster_waning <- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_reinf        <- array(dim = c(num_runs, num_days_sim, 11))
  mean_susceptible_vac_reinfvac <- array(dim = c(num_runs, num_days_sim, 11))
  mean_vaccinated_age  <- array(dim = c(num_runs, num_days_sim, 11))

  scen_hosp_adm     <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  scen_hosp_adforcovid     <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  scen_Rt_infection <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_Rt_mild_sev  <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_Rt_hosp      <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_new_infect   <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_new_reinfect <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_new_sympt    <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_mortality    <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_VOC_mild_alpha     <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_VOC_mild_delta     <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_VOC_mild_omicron   <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  scen_VOC_mild_ba4ba5    <- matrix(NA, nrow= num_days_sim, ncol =num_runs)
  
  scen_hosp_load    <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  scen_ICU_load     <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  scen_ICU_adm      <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 

  scen_hosp_exit    <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
    
  nvac_hosp_adm      <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  nvac_age_hosp_adm  <- array(dim = c(num_runs, num_days_sim, 11))
  new_hosp_icu_vac_d2 <- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  new_hosp_icu_vac_age_d2 <- array(dim = c(num_runs, num_days_sim, 11))
  new_hosp_icu_vac_booster<- matrix(NA, nrow= num_days_sim, ncol =num_runs) 
  
  mean_rna1_vaccinated    <- array(dim = c(num_runs, num_days_sim, 11)) 
  mean_rna2_vaccinated    <- array(dim = c(num_runs, num_days_sim, 11)) 
  mean_adeno1_vaccinated  <- array(dim = c(num_runs, num_days_sim, 11))
  mean_adeno2_vaccinated  <- array(dim = c(num_runs, num_days_sim, 11))
  mean_waning_vaccinated  <- array(dim = c(num_runs, num_days_sim, 11))
  
  i <- 1
  time_stamp_loop <- Sys.time()
  for (i in 1:num_runs){
    set.seed(20210101 + i)
    
    smd_print_progress(i,num_runs,time_stamp_loop=time_stamp_loop,par_nodes_info = par_nodes_info)
    run_param <- unlist(parms_chains[nrow(parms_chains)-(i_seq[i]-1), ])
    
    pred_all <- run_model(parms = run_param, 
                       contact_data = contact_data, 
                       method = "stochastic", 
                       ndays_sim = num_days_sim,
                       vaccine_uptake = vacc_schedule,
                       be_ref_data = be_ref_data,
                       db_cnt_adjust = db_cnt_adjust)
    
    # hospital admissions
    pred_hosp1          <- pred_all$total_new_hosp_icu
    hosp_age_adm[i,,]   <- as.matrix(pred_hosp1)       
    scen_hosp_adm[,i]   <- as.matrix(rowSums(pred_all$total_new_hosp_icu[,-1]))
    
    # hospital admissions for covid (with estimation/correction)
    pred_hosp1_for_covid       <- pred_all$total_new_hosp_icu_for_covid
    hosp_age_adforcovid[i,,]   <- as.matrix(pred_hosp1_for_covid)       
    scen_hosp_adforcovid[,i]   <- as.matrix(rowSums(pred_all$total_new_hosp_icu_for_covid[,-1]))
  
    # infections
    new_cases           <- pred_all$total_new_infections
    scen_cases[i,,]     <- as.matrix(new_cases)      
    scen_new_infect[,i] <- as.matrix(rowSums(new_cases[,-1])) # new_cases
    
    # reinfections
    new_reinfect          <- pred_all$total_new_reinfections
    scen_new_reinfect[,i] <- as.matrix(rowSums(new_reinfect[,-1])) # new reinfections
    
    # mild cases
    new_mild_cases        <- pred_all$total_new_mild_infections
    scen_mild_cases[i,,]  <- as.matrix(new_mild_cases)  
    scen_new_sympt[,i]    <- as.matrix(rowSums(new_mild_cases[,-1])) # new_cases
    
    # hospital and icu load 
    scen_hosp_load[,i]   <- as.matrix(rowSums(pred_all$hosp_icu_load[,-1]))
    scen_ICU_load[,i]    <- as.matrix(rowSums(pred_all$icu_load[,-1]))
    scen_ICU_adm[,i]     <- as.matrix(rowSums(pred_all$total_new_icu[,-1]))

    # prevalence and incidence (DALY)
    prev_E_sims[i,,]         <- as.matrix(pred_all$prev_E[,-1])
    prev_I_presymp_sims[i,,] <- as.matrix(pred_all$prev_I_presym[,-1])
    prev_I_asymp_sims[i,,]   <- as.matrix(pred_all$prev_I_asym[,-1])
    prev_I_mild_sims[i,,]    <- as.matrix(pred_all$prev_I_mild[,-1])
    prev_I_severe_sims[i,,]  <- as.matrix(pred_all$prev_I_sev[,-1])
    prev_I_hosp_sims[i,,]    <- as.matrix(pred_all$hosp_icu_load[,-1])
    prev_I_icu_sims[i,,]     <- as.matrix(pred_all$icu_load[,-1])
    prev_M_sims[i,,]         <- as.matrix(apply(pred_all$total_new_deaths[,-1],2,cumsum))
    prev_R_sims[i,,]         <- as.matrix(pred_all$prev_R[,-1])
    
    inc_I_hosp_sims[i,,]     <- as.matrix(pred_all$total_new_hosp_icu[,-1])
    inc_I_icu_sims[i,,]      <- as.matrix(pred_all$total_new_icu[,-1])
    
    inc_I_hosp_for_covid_sims[i,,]     <- as.matrix(pred_all$total_new_hosp_icu_for_covid[,-1])
    # hospital exit
    scen_hosp_exit[,i]   <- as.matrix(rowSums(pred_all$total_hosp_exit[,-1]))
    
    # Rt
    scen_Rt_infection[,i] <- get_Rt(vect_cases = rowSums(new_cases[,-1]),
                                    vect_dates = new_cases$day)$Rt
    scen_Rt_mild_sev[,i] <- get_Rt(vect_cases = rowSums(new_mild_cases[,-1]),
                                  vect_dates = new_mild_cases$day)$Rt
    scen_Rt_hosp[,i]     <- get_Rt(vect_cases = rowSums(pred_hosp1[,-1]),
                                   vect_dates = pred_hosp1[,1])$Rt
    
    # mortality
    scen_mortality[,i]   <- as.matrix(rowSums(pred_all$total_new_deaths[,-1]))
    mortality_age[i,,]   <- as.matrix(pred_all$total_new_deaths) 
    
    # make sure that "VOC_ba4ba5_start" is part of run_param
    if(!"VOC_ba4ba5_start" %in% names(run_param)){
      run_param['VOC_ba4ba5_start'] <- num_days_sim
    }
    
    # VOC
    d_alpha   <- 1:min(run_param['VOC_omicron_start'],num_days_sim)
    d_delta   <- min(run_param['VOC_delta_start'],num_days_sim):min(run_param['VOC_ba4ba5_start'],num_days_sim)
    d_omicron <- min(run_param['VOC_omicron_start'],num_days_sim):num_days_sim
    d_ba4ba5  <- min(run_param['VOC_ba4ba5_start'],num_days_sim):num_days_sim
    scen_VOC_mild_alpha[d_alpha,i]     <- as.matrix(pred_all$p_VOC)[d_alpha]
    scen_VOC_mild_delta[d_delta,i]     <- as.matrix(1-pred_all$p_VOC)[d_delta]
    scen_VOC_mild_omicron[d_omicron,i] <- as.matrix(pred_all$p_VOC)[d_omicron]
    scen_VOC_mild_ba4ba5[d_ba4ba5,i]   <- as.matrix(1-pred_all$p_VOC)[d_ba4ba5]
    
    # cumulative infections
    scen_cases_cum     <- apply(as.matrix(new_cases),2,cumsum)
    scen_cases_cum[,1] <- as.matrix(new_cases)[,1]
    cumul_cases[i,,]   <- scen_cases_cum
    
    # susceptible
    mean_susceptible_full[i,,]          <- as.matrix(pred_all$mean_susceptible_full)  
    mean_susceptible_vac_d1[i,,]        <- as.matrix(pred_all$mean_susceptible_vac_d1)  
    mean_susceptible_vac_d2[i,,]        <- as.matrix(pred_all$mean_susceptible_vac_d2) 
    mean_susceptible_vac_d2_waning[i,,] <- as.matrix(pred_all$mean_susceptible_vac_d2_waning)  
    mean_susceptible_vac_booster[i,,]   <- as.matrix(pred_all$mean_susceptible_vac_booster)  
    mean_susceptible_vac_booster_waning[i,,]   <- as.matrix(pred_all$mean_susceptible_vac_booster_waning)  
    mean_susceptible_reinf[i,,]         <- as.matrix(pred_all$mean_susceptible_reinf)  
    mean_susceptible_vac_reinfvac[i,,]  <- as.matrix(pred_all$mean_susceptible_vac_reinfvac)  
    
    mean_vaccinated_age[i,,]       <- as.matrix(pred_all$mean_vaccinated_age)
  
    # non-vaccinated hosp admissions
    nvac_hosp_adm[,i]   <- as.matrix(rowSums(pred_all$total_new_hosp_icu_nvac[,-1]))
    nvac_age_hosp_adm[i,,]   <- as.matrix(pred_all$total_new_hosp_icu_nvac)
    new_hosp_icu_vac_d2[,i] <- as.matrix(rowSums(pred_all$total_new_hosp_icu_vac_d2[,-1]))
    new_hosp_icu_vac_age_d2[i,,] <- as.matrix(pred_all$total_new_hosp_icu_vac_d2)
    new_hosp_icu_vac_booster[,i] <- as.matrix(rowSums(pred_all$total_new_hosp_icu_vac_booster[,-1]))
    
    # vaccinated
    flag_rna     <- grepl('rna',names(pred_all$mean_dose_vaccinated))
    flag_adeno   <- grepl('adeno',names(pred_all$mean_dose_vaccinated))  
    flag_d1      <- grepl('d1',names(pred_all$mean_dose_vaccinated))
    flag_d2      <- grepl('d2',names(pred_all$mean_dose_vaccinated))
    flag_booster <- grepl('booster',names(pred_all$mean_dose_vaccinated)) & ! grepl('waning',names(pred_all$mean_dose_vaccinated))
    
    mean_rna1_vaccinated[i,,]    <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_rna & flag_d1))])
    mean_rna2_vaccinated[i,,]    <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_rna & flag_d2))])
    mean_adeno1_vaccinated[i,,]  <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_adeno & flag_d1))])
    mean_adeno2_vaccinated[i,,]  <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_adeno & flag_d2))])
    mean_waning_vaccinated[i,,]  <- as.matrix(pred_all$mean_dose_vaccinated[,c(1,which(flag_booster))])
  }
  
  output_file_name_generic <- paste0(output_dir,'XXX',gsub('\\.','_incr_p',cnt_adjust_tag),'_',vacc_schedule_tag,'.rds')
  
  saveRDS(hosp_age_adm,       file=gsub('XXX','scenario1_incr',output_file_name_generic))
  saveRDS(hosp_age_adm,       file=gsub('XXX','hosp_age_adm_incr',output_file_name_generic))
  saveRDS(hosp_age_adforcovid,file=gsub('XXX','hosp_age_adforcovid_incr',output_file_name_generic))
  saveRDS(scen_cases,         file=gsub('XXX','scen_cases_incr',output_file_name_generic))
  saveRDS(scen_mild_cases,    file=gsub('XXX','scen_mild_cases_incr',output_file_name_generic))
  saveRDS(cumul_cases,        file=gsub('XXX','cumul_cases_incr',output_file_name_generic))
  saveRDS(mean_susceptible_full,              file=gsub('XXX','mean_susceptible_full_incr',output_file_name_generic))
  saveRDS(mean_susceptible_vac_d1,            file=gsub('XXX','mean_susceptible_vac_d1_incr',output_file_name_generic))
  saveRDS(mean_susceptible_vac_d2,            file=gsub('XXX','mean_susceptible_vac_d2_incr',output_file_name_generic))
  saveRDS(mean_susceptible_vac_d2_waning,     file=gsub('XXX','mean_susceptible_vac_d2_waning_incr',output_file_name_generic))
  saveRDS(mean_susceptible_vac_booster,       file=gsub('XXX','mean_susceptible_vac_booster_incr',output_file_name_generic))
  saveRDS(mean_susceptible_vac_booster_waning,file=gsub('XXX','mean_susceptible_vac_booster_waning_incr',output_file_name_generic))
  saveRDS(mean_susceptible_reinf,             file=gsub('XXX','mean_susceptible_reinf_incr',output_file_name_generic))
  saveRDS(mean_susceptible_vac_reinfvac,      file=gsub('XXX','mean_susceptible_vac_reinfvac_incr',output_file_name_generic))
  saveRDS(mean_vaccinated_age,                file=gsub('XXX','mean_vaccinated_age_incr',output_file_name_generic))
  
  saveRDS(mean_rna1_vaccinated,  file=gsub('XXX','mean_rna1_vaccinated_incr',output_file_name_generic))
  saveRDS(mean_rna2_vaccinated,  file=gsub('XXX','mean_rna2_vaccinated_incr',output_file_name_generic))
  saveRDS(mean_adeno1_vaccinated,file=gsub('XXX','mean_adeno1_vaccinated_incr',output_file_name_generic))
  saveRDS(mean_adeno2_vaccinated,file=gsub('XXX','mean_adeno2_vaccinated_incr',output_file_name_generic))
  saveRDS(mean_waning_vaccinated,file=gsub('XXX','mean_waning_vaccinated_incr',output_file_name_generic))
  
  saveRDS(scen_new_infect,      file=gsub('XXX','new_infect',output_file_name_generic))
  saveRDS(scen_new_reinfect,    file=gsub('XXX','new_reinfect',output_file_name_generic))
  saveRDS(scen_new_sympt,       file=gsub('XXX','new_sympt',output_file_name_generic))
  saveRDS(scen_hosp_load,       file=gsub('XXX','hosp_load_incr',output_file_name_generic))
  saveRDS(scen_Rt_infection,    file=gsub('XXX','scen_Rt_infection_incr',output_file_name_generic))
  saveRDS(scen_Rt_mild_sev,     file=gsub('XXX','scen_Rt_mild_sev_incr',output_file_name_generic))
  saveRDS(scen_Rt_hosp,         file=gsub('XXX','scen_Rt_hosp_incr',output_file_name_generic))
  saveRDS(scen_hosp_adm,        file=gsub('XXX','hosp_adm_incr',output_file_name_generic))
  saveRDS(scen_hosp_adforcovid, file=gsub('XXX','hosp_adforcovid_incr',output_file_name_generic))
  saveRDS(scen_ICU_load,        file=gsub('XXX','icu_load_incr',output_file_name_generic))
  saveRDS(scen_ICU_adm,         file=gsub('XXX','icu_adm_incr',output_file_name_generic))
  saveRDS(scen_hosp_exit,       file=gsub('XXX','hosp_exit_incr',output_file_name_generic))
  saveRDS(scen_mortality,       file=gsub('XXX','scen_mortality',output_file_name_generic))
  saveRDS(mortality_age,        file=gsub('XXX','mortality_age_incr',output_file_name_generic))
  saveRDS(scen_VOC_mild_alpha,  file=gsub('XXX','scen_VOC_mild_alpha_incr',output_file_name_generic))
  saveRDS(scen_VOC_mild_delta,  file=gsub('XXX','scen_VOC_mild_delta_incr',output_file_name_generic))
  saveRDS(scen_VOC_mild_omicron,file=gsub('XXX','scen_VOC_mild_omicron_incr',output_file_name_generic))
  saveRDS(scen_VOC_mild_ba4ba5, file=gsub('XXX','scen_VOC_mild_ba4ba5_incr',output_file_name_generic))

  # burden of disease figures
  saveRDS(prev_E_sims,         file=gsub('XXX','prev_E_sims_incr',output_file_name_generic))
  saveRDS(prev_I_presymp_sims, file=gsub('XXX','prev_I_presymp_sims_incr',output_file_name_generic))
  saveRDS(prev_I_asymp_sims,   file=gsub('XXX','prev_I_asymp_sims_incr',output_file_name_generic))
  saveRDS(prev_I_mild_sims,    file=gsub('XXX','prev_I_mild_sims_incr',output_file_name_generic))
  saveRDS(prev_I_severe_sims,  file=gsub('XXX','prev_I_severe_sims_incr',output_file_name_generic))
  saveRDS(prev_I_hosp_sims,    file=gsub('XXX','prev_I_hosp_sims_incr',output_file_name_generic))
  saveRDS(prev_R_sims,         file=gsub('XXX','prev_R_sims_incr',output_file_name_generic))
  saveRDS(prev_M_sims,         file=gsub('XXX','prev_M_sims_incr',output_file_name_generic))
  saveRDS(prev_I_icu_sims,     file=gsub('XXX','prev_I_icu_sims_incr',output_file_name_generic))
  saveRDS(inc_I_hosp_sims,     file=gsub('XXX','inc_I_hosp_sims_incr',output_file_name_generic))
  saveRDS(inc_I_icu_sims,      file=gsub('XXX','inc_I_icu_sims_incr',output_file_name_generic))
  saveRDS(inc_I_hosp_for_covid_sims,file=gsub('XXX','inc_I_hosp_for_covid_sims_incr',output_file_name_generic))
  
  vac_hosp_adm <- scen_hosp_adm - nvac_hosp_adm
  saveRDS(nvac_hosp_adm, file=gsub('XXX','nvac_hosp_adm_incr',output_file_name_generic))
  saveRDS(vac_hosp_adm,  file=gsub('XXX','vac_hosp_adm_incr',output_file_name_generic))
  
  vac_age_hosp_adm <- hosp_age_adm - nvac_age_hosp_adm
  vac_age_hosp_adm[,,1] <- hosp_age_adm[,,1]
  saveRDS(vac_age_hosp_adm,  file=gsub('XXX','vac_age_hosp_adm_incr',output_file_name_generic))
  saveRDS(nvac_age_hosp_adm, file=gsub('XXX','nvac_age_hosp_adm_incr',output_file_name_generic))
  
  saveRDS(new_hosp_icu_vac_d2,      file=gsub('XXX','new_hosp_adm_vac_d2_incr',output_file_name_generic))
  saveRDS(new_hosp_icu_vac_age_d2,  file=gsub('XXX','new_hosp_adm_vac_age_d2_incr',output_file_name_generic))
  saveRDS(new_hosp_icu_vac_booster, file=gsub('XXX','new_hosp_adm_vac_booster_incr',output_file_name_generic))

  if(!all(is.na(vacc_schedule))){
    # save vaccine uptake (given the number of days simulated)
    # note: this is not the same as the 'V_mat' in the model kernel, which represents vaccine protection
    uptake_doses_total <- colSums(vacc_schedule[vacc_schedule$date < (get_scm_start_date()+num_days_sim),-1])
    saveRDS(uptake_doses_total,file=gsub('XXX','uptake_doses_total_incr',output_file_name_generic))
  }
  
  ## plot results ----
  smd_print("Save output files")
  pdf(file=gsub('rds','pdf',gsub('XXX','scenario_cases_mean_incr',output_file_name_generic)),9,6)
  par(mar=c(4.5,5,1,1))
  scm_scenario_cp <- as.Date(db_cnt_adjust$date)[db_cnt_adjust$value>0] 
  x_axis_scm_day = c(0,num_days_sim) 
  x_axis_vaccine = c(0,num_days_sim)
  
  output_files <- dir(output_dir,full.names = T,pattern = paste(vacc_schedule_tag))
  multi_plot_incidence_time(output_files       = output_files,
                            bool_polygon_multi = T,
                            col_lib            = col_lib,
                            be_ref_data        = be_ref_data,
                            db_C_sim           = contact_data$db_C_sim,
                            default_y_axis     = default_y_axis,
                            x_axis_scm_day     = x_axis_scm_day,
                            scm_scenario_cp    = scm_scenario_cp)
  dev.off()
  
  # print 1 plot
  smd_print("Return additional hospital admission plot")
  multi_plot_incidence_time(output_files       = output_files,
                            output_tag_multi   ='/hosp_adm',
                            bool_polygon_multi = T,
                            col_lib            = col_lib,
                            be_ref_data        = be_ref_data,
                            db_C_sim           = contact_data$db_C_sim,
                            default_y_axis     = default_y_axis,
                            x_axis_scm_day     = x_axis_scm_day,
                            scm_scenario_cp    = scm_scenario_cp)
  
  # # optional: plot age-specific results
  # pdf(file=gsub('rds','pdf',gsub('XXX','scenario_cases_mean_age',output_file_name_generic)),9,6)
  # multi_plot_incidence_age_time(output_files   = output_files,
  #                               col_lib        = col_lib,
  #                               x_axis_scm_day  = x_axis_scm_day,
  #                               x_axis_vaccine = x_axis_vaccine)
  # 
  # dev.off()
  
} # end parallel foreach

smd_stop_cluster()
}
