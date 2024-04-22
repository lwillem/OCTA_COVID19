########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains many helper function to run the stochastic model, 
# to explore the results and estimate model parameters.
#
# Copyright 2024, SIMID                                        
########################################################################### #

############################################################################ #
# PROVIDE PUBLIC HEALTH AGENCY DATA  ----                                 
############################################################################ #

select_regional_data <- function(ref_data,sel_region){
  
  if(is.na(sel_region) | sel_region == 'belgium'){
    return(ref_data)
  }
  
  if(!any(names(ref_data) == 'REGION')){
    print("NO REGION SPECIFIC REFERENCE DATA")
    return(ref_data)
  }
  
  # else
  return(ref_data[!is.na(ref_data$REGION) & grepl(substr(sel_region,1,5),tolower(ref_data$REGION)),])
}


get_be_pop_matrix <- function(n_row){
  
  pop_be2020_all <- read.table('data/pop_be2020_statbel.csv',sep=',',header=T)
  pop_be2020 <- pop_be2020_all[,3]
  
  return(matrix(rep(pop_be2020, each = n_row),
                nrow = n_row,
                ncol = length(pop_be2020)))
}

get_regional_pop <- function(region = 'belgium',bool_10y = TRUE){
  
  pop_be2020_all <- read.table('data/pop_be2020_statbel_download.csv',sep=',',header=T)
  unique(pop_be2020_all$Region)
  if(is.na(region) | region == 'belgium'){
    pop_be2020_all <- pop_be2020_all[nchar(pop_be2020_all$Region)==0,]
  } else{
    pop_be2020_all <- pop_be2020_all[grepl(substr(region,1,5),tolower(pop_be2020_all$Region)),]
  }
  
  if(length(unique(pop_be2020_all$Region))>1){
    warning('provided "region" is not vallid, use belgium, flanders, walloon or brussels')
    return(NA)
  }
  
  # aggregate 90-99 and +100
  pop_be2020_all$Population.on.January.1st.2020[nrow(pop_be2020_all)-1] <- sum(pop_be2020_all$Population.on.January.1st.2020[nrow(pop_be2020_all)-(1:0)])
  pop_be2020_all <- pop_be2020_all[-nrow(pop_be2020_all),]
  
  # optional: aggregate 5y age groups into 10y age groups
  if(bool_10y){
    pop2020 <- colSums(matrix(c(pop_be2020_all$Population.on.January.1st.2020),ncol=10,nrow=2,byrow=F))
    return(pop2020)
  } else{
    return(pop_be2020_all)
  }

}

get_latest_incidence_data <- function(sel_region='belgium'){
  
  if(!sel_region %in% get_region()){
    smd_print(c('The provided region ',sel_region,' is not part of the specified regions:', get_region()),WARNING = TRUE)
  }
  
  # get file name of stored reference data
  ref_data_file_name <- dir('data',pattern=paste0('covid19_reference_data_.*',sel_region,'.*csv'),full.names = T)
  
  # if multiple exists, use the most recent
  ref_data_file_name <- ref_data_file_name[length(ref_data_file_name)]
  
  ref_data      <- read.table(ref_data_file_name,sep=',',header=T)
  ref_data$date <- as.Date(ref_data$date)
  
  return(ref_data)
}


# Function to efficiently generate summary tables
get_age_table <- function(data_transm){
  
  # overal summary
  summary_table_general         <- aggregate(formula('CASES ~ DATE'), data=data_transm, sum)
  names(summary_table_general)  <- c('date','cases')
  
  # specific summary
  summary_table                  <- dcast(data_transm, formula(paste('DATE', '~' ,'AGEGROUP')), value.var='CASES', sum)
  
  # update names
  names(summary_table)           <- c('date',paste('cases',names(summary_table)[-1],sep='_'))
  
  # merge
  summary_table <- merge(summary_table,summary_table_general)
  
  # omit NA's (default on LW's MACOS but not on VSC cluster)
  # R versions: 3.5.3 (MACOS) vs. 3.5.1 (VSC)
  summary_table <- na.omit(summary_table,cols='date')
  
  #check
  head(summary_table)
  tail(summary_table)
  
  # return
  return(summary_table)
}

get_Rt <- function(vect_cases,vect_dates){
  
  # check for NA
  flag_na <- is.na(vect_cases)
  
  if(any(vect_cases<0)){
    warning('Negative number of cases! Not possible to calculate Rt')
    return(list(Rt=NA))
  }
  
  # use all data, if the last observation needs to be removed, this should be done in previous data cleaning
  cont.series = vect_cases
  cont.dates  = vect_dates
  
  # join
  dta.cont=data.frame(cont.dates,cont.series)
  names(dta.cont) <- c("dates", "I")
  
  # make sure the dates are consecutive
  dates.dummy <- data.frame(dates=seq(min(dta.cont$dates),max(dta.cont$dates),1))
  dta.cont  <- merge(dates.dummy,dta.cont,all.x = T)
  dta.cont$I[is.na(dta.cont$I)] <- 0 # set NA to "0"
  
  ## calculate effective reproduction number
  suppressMessages(
    Rt_mod <-
        estimate_R(
          dta.cont,
          method = "parametric_si",
          config = make_config(list(mean_si = 4.7, 
                                    std_si = 2.9)))
   )
  
    # add the reproduction number over the 7-day window finishing on that day.
    dta.cont$Rt <- NA
    dta.cont$Rt[Rt_mod$R$t_end] <- Rt_mod$R$`Mean(R)`
  
    # account for last observation
    dta.cont[nrow(dta.cont)+1,] <- dta.cont[nrow(dta.cont),]+1
    dta.cont[nrow(dta.cont),-1] <- NA
  
    # remove imputed dates
    if(length(vect_dates) < nrow(dta.cont)){
      dta.cont <- dta.cont[dta.cont$dates %in% vect_dates,]
    }
    
  return(dta.cont)
}

get_Rt_epilps <- function(vect_cases,vect_dates){
  
  # defensive programming
  if(!"EpiLPS" %in% installed.packages()){
    return(NA)
  }
  
  # check for NA
  flag_na <- is.na(vect_cases)
  
  if(any(vect_cases<0)){
    warning('Negative number of cases! Not possible to calculate Rt')
    return(list(Rt=NA))
  }
  
  # use all data, if the last observation needs to be removed, this should be done in previous data cleaning
  cont.series = vect_cases
  cont.dates  = vect_dates
  
  # mean 3 days and sd 2.9
  sidistr=discr_si(k=seq(0, 14), mu=3, sigma=2.9);sidistr=sidistr/sum(sidistr)
  
  # calculation Rt
  dta.rt <- epilps(incidence = cont.series, 
                     serial_interval = sidistr, verbose = FALSE)
  
  return(dta.rt)
}

recap_Rt <- function(output_dir,tag_incidence,tag_output){
  
  # get file name
  rds_filenames <- dir(output_dir,recursive = T,full.names = T)
  rds_filenames <- rds_filename_incidence[grepl(tag_incidence,rds_filename_incidence)]
  
  for(rds_filename_incidence in rds_filenames){
    # load data
    incidence <- readRDS(rds_filename_incidence)
    Rt_new    <- incidence*NA
    dim(incidence)
    for(i_sim in 1:ncol(incidence)){
      Rt_new[,i_sim] <- get_Rt(incidence[,i_sim],
                               1:nrow(incidence))$Rt
    }  
    
    # get output file name
    rds_filename_Rt <- gsub(tag_incidence,tag_output,rds_filename_incidence)
    
    # save output
    saveRDS(Rt_new,file = rds_filename_Rt)
    print(rds_filename_Rt)
  }
}

## calculate exponential growth rate
get_growth_rate <- function(f_obs, exp_window = 7){
  
  # exp_window <- 7  # days to calculate the growth
  out_r    = rep(NA,length(f_obs))
  
  for (i in exp_window:length(f_obs)){
    temp_df = as.data.frame(cbind(1:exp_window,f_obs[(i-exp_window+1):i]))
    names(temp_df)=c("day","cases")
    mylm = lm(log10(cases+1) ~ day, temp_df)
    out_r[i] <- coefficients(mylm)[[2]]
  }
  
  return(out_r)
}


get_age_groups <- function(){
  age_breaks      <- seq(1,110,10)
  age_start       <- age_breaks[-length(age_breaks)]
  age_end         <- age_breaks[-1]-1
  age_tags        <- paste(age_start-1,age_end-1,sep="_")
  return(age_tags)
}

get_age_group_labels <- function(){
  return(gsub('_','-',get_age_groups()))
}

# CONVERT A CALENDAR DATE INTO A SIMULATION DAY INDEX
# default date format: YYYY-MM-DD
#note: March 1st, 2020 == "day 0"
sim_date2day <- function(x_date,x_format='%Y-%m-%d'){
  return(as.numeric(as.Date(as.character(x_date),x_format) - as.Date('2020-03-01')))
}

# CONVERT SIMULATION DAY INDEX INTO CALENDAR DATE
#note: March 1st, 2020 == "day 0"
sim_day2date <- function(x_day_index){
  return(as.Date('2020-03-01') + x_day_index)
}

# CONVERT SIMULATION TIME STEP INTO DAY INDEX
#note: the 1st step in 2020 == "step 1"
sim_step2day <- function(time_step){
  # return(floor(time_step/24))
  return(floor((time_step-1)/4))
}

# CONVERT SIMULATION TIME STEP INTO CALENDAR DATE
sim_step2date <- function(time_step){
  return(sim_day2date(sim_step2day(time_step)))
}


## Helper functions ----
##----------------- -

add_date_axis <- function(f_dates, side = 1, num_intervals = 12, bool_grid=FALSE){
  
  x_ticks <- pretty(as.Date(f_dates),num_intervals)
  x_labels <- format(x_ticks,'%d/%m')
  axis(side,x_ticks,x_labels)
  abline(v=x_ticks,lty=3,col='grey')
}
  
# scenario_data <- scenarioXa; col <- "blue"
add_polygon <- function(scenario_data,scenario_dates,col,alpha_val = 0.5,mean_lty=1){
  
  # make sure that "scenario_data" is a data.frame
  if(typeof(scenario_data)=="double"){
    scenario_data <- data.frame(scenario_data)
  }
  
  # get column id's
  id_mean  <- grepl('_mean',names(scenario_data))
  id_lower <- grepl('_LL',names(scenario_data))
  id_upper <- grepl('_UL',names(scenario_data))
  
  id_975  <- grepl('_975',names(scenario_data))
  id_025  <- grepl('_025',names(scenario_data))
  
  # remove NA's
  flag_na <- is.na(rowSums(scenario_data))
  if(any(flag_na)){
    scenario_data <- scenario_data[!flag_na,]
    scenario_dates <- scenario_dates[!flag_na]
  }
  
  # make sure that the dates are sorted
  scenario_data  <- scenario_data[order(scenario_dates),]
  scenario_dates <- scenario_dates[order(scenario_dates)]
  
  # add uncertainty interval
  polygon(x = c(scenario_dates,rev(scenario_dates)),
          y = c(scenario_data[,id_lower],rev(scenario_data[,id_upper])),
          col = alpha(col,alpha_val),
          border = NA)
  
  # add mean
  lines(x = scenario_dates,
        y = scenario_data[,id_mean],
        lty = mean_lty,
        col = col,
        lwd = 3)
  
  if(any(id_975) && any(id_025)){
    lines(x = scenario_dates,
          y = scenario_data[,id_975],
          lty = 2,
          col = col,
          lwd = 1)
    lines(x = scenario_dates,
          y = scenario_data[,id_025],
          lty = 2,
          col = col,
          lwd = 1)
  }
  
}

# add vertical line on given date + label on x-axis
add_vertical_line <- function(date_string,bool_text,date_tag = '',pos_factor = 0.05){
  
  plot_limits <- par("usr")
  
  v_date <- as.Date(date_string)
  abline(v=v_date,lty=2)
  #axis(1,v_date,format(v_date,'%d/%m'),
  #cex.axis=0.5,padj=-3,tck=-0.005)
  
  if(bool_text)
  {
    v_text <- ifelse(nchar(date_tag)>0,date_tag,format(v_date,'%d/%m'))
    text(x = v_date-1,
         y = mean(plot_limits[3:4])*pos_factor,
         #paste(format(v_date,'%d/%m'),date_tag),
         v_text,
         srt=90, pos=3, offset = +1.5,cex=0.6)
  }
}


# add copyright statement to the bottom right corner of your figure.
add_copyright <- function(text.cex = 0.4,bool_left = FALSE){
  plot_limits <- par("usr")
  x_value <- ifelse(bool_left,plot_limits[1],plot_limits[2])
  y_value <- plot_limits[3] + diff(plot_limits[3:4])/50
  x_adj <- ifelse(bool_left,-0.1,1)
  text(x = x_value,
       y = y_value,
       labels = ('© SIMID / UAntwerpen / UHasselt  '),
       cex = text.cex,
       adj = x_adj)
}

# return model version as included in the README file 
get_scm_version <- function(filename_readme = './README.md'){
  
  readme <- readLines(filename_readme)
  str_version <- readme[grep("Version\\:", readme)]
  str_version <- unlist(strsplit(str_version,' '))
  
  return(as.double(str_version[length(str_version)]))
}

# set SCM start date
get_scm_start_date <- function(){
  return(as.Date('2020-03-01'))
}

# parm_name <- i_param
get_ndays_sim <- function(parm_name,ndays_sim,extra_days=40){
  
  if(is.na(parm_name)){
    return(ndays_sim)
  }
  
  pname <- names(parm_name)
  if(!grepl('coef_w',pname)){
    return(ndays_sim)
  }
  wave_id <- gsub('.*coef_w','',pname)
  wave_id <- gsub('_.*','',wave_id)
  
  cp_wave <- ifelse(grepl('comix_coef_w',pname),
                    get_CoMix_change_day(as.numeric(wave_id)),
                    get_add_change_day(as.numeric(wave_id)))
  
  return(min(ndays_sim,cp_wave+extra_days))
}

adjust_q_param <- function(parms,parms_names,parms_names_estim,
                           bool_setup = FALSE, 
                           opt_waves = 9:30){
  
  if(!bool_setup && !any(grepl('stage_aggr',parms_names_estim))){
    return(parms)
  }
  
  if(typeof(opt_waves) != 'list'){
    tbl_q_waves <- list(opt_waves)
  } else{
    tbl_q_waves <- opt_waves
  }
  
  length(tbl_q_waves)
  
  # option: initial set-up of aggregated q-param
  if(bool_setup){
    for(i_q in 1:length(tbl_q_waves)){
      q_aggr_tag <- paste0('stage_aggr',i_q,'_age')
      sel_wave_init  <- tbl_q_waves[[i_q]][1]
      parms[paste0(q_aggr_tag,1:10)] <- parms[get_colnames(parms_names,tag_list = paste0('stage',sel_wave_init,'_age'))]
    }
  }
  
  # copy aggregated q-parameters to wave-specific q-parameters
  for(i_q in 1:length(tbl_q_waves)){
    q_aggr_tag <- paste0('stage_aggr',i_q,'_age')
    if(any(grepl(q_aggr_tag,parms_names_estim))){
      sel_waves <- tbl_q_waves[[i_q]]
      parms[get_colnames(parms_names,tag_list = paste0('stage',sel_waves,'_'))]  <- rep(parms[paste0(q_aggr_tag,1:10)],length(sel_waves))
    } 
  }
  
  return(parms)
}


## Region specific functions   ----
##---------------------- -

db_region <- c("belgium","brussels","flanders","wallonia")
get_region_id <- function(str_region){
  return(which(db_region == str_region))
}

get_region <- function(id_region){
  return(db_region[id_region])
}


## Uptake functions   ----
##---------------------- -

# note: time_step is best expressed in "days" (or at least the same unit)
approx_uptake <- function(aggr_uptake, target_time_step, dose_col=NA,max_date = NA){
  
  # set output dates
  x_out <- seq(min(aggr_uptake$date),max(c(aggr_uptake$date,max_date),na.rm = T),target_time_step)
  
  # get current aggregation
  aggr_time_steps <- as.numeric(c(diff(aggr_uptake$date),1))
  
  # prepare final data.frame and add dates
  if(!any(is.na(dose_col))){
    target_uptake <- data.frame(matrix(0,nrow=length(x_out),ncol=(1+length(dose_col)*2)))
    names(target_uptake) <- c('date',
                              paste0('adeno_',dose_col),
                              paste0('rna_',dose_col))
    target_uptake$date   <- x_out
  } else{
    target_uptake        <- data.frame(matrix(0,nrow=length(x_out),ncol=ncol(aggr_uptake)))
    names(target_uptake) <- names(aggr_uptake)
    target_uptake$date   <- x_out
  }

  sel_column <- names(aggr_uptake)[2]
  for(sel_column in names(aggr_uptake)[-1]){
    
    aggr_uptake[,sel_column]
    aggr_uptake$date
    
    # linear extrapolate the doses over the days of the month
    approx(x = aggr_uptake$date,
           y = aggr_uptake[,sel_column] / aggr_time_steps * target_time_step,
           xout = x_out,
           yright = aggr_uptake[nrow(aggr_uptake),sel_column] / aggr_time_steps[length(aggr_time_steps)] * target_time_step,
           method='constant') -> aggr_uptake_hour
    
    target_uptake[,sel_column]   <- aggr_uptake_hour$y
  }
  
  return(target_uptake)
}

uptake2protection <- function(uptake_time,
                              delay_protection_rna,
                              delay_protection_adeno){
  
  # format and make sure to work with unique date ids
  uptake_time$date        <- as.Date(uptake_time$date)

  if(length(unique(uptake_time$date)) != nrow(uptake_time)){
    # (re)set output dates (ERROR PRONE!)
    target_time_step        <- sum(uptake_time$date[1] == uptake_time$date)
    uptake_time$date <- seq(min(uptake_time$date),max(uptake_time$date),1/target_time_step)
  }
  
  # select uptake: rna-based vaccine
  uptake_time_rna      <- uptake_time[,grepl('rna',names(uptake_time))]
  uptake_time_rna$date <- uptake_time$date + delay_protection_rna

  # select uptake: adeno-based vaccine
  uptake_time_adeno      <- uptake_time[,grepl('adeno',names(uptake_time))]
  uptake_time_adeno$date <- uptake_time$date + delay_protection_adeno

  # merge datasets
  uptake_time_delay <- data.frame(date       = c(uptake_time$date[-nrow(uptake_time)],
                                                 seq(max(uptake_time$date),max(uptake_time$date)+30,1/4)),
                                  dummy = 1)
  
  # use numeric 'step' instead of 'date' to prevent merge issues
  uptake_time_rna$date_num   <- as.numeric(uptake_time_rna$date)
  uptake_time_adeno$date_num <- as.numeric(uptake_time_adeno$date)
  uptake_time_delay$date_num <- as.numeric(uptake_time_delay$date)
  uptake_time_rna$date   <- NULL
  uptake_time_adeno$date <- NULL
  
  # merge
  uptake_time_delay <- merge(uptake_time_delay,uptake_time_adeno,by='date_num',all.x = T,all.y=F)
  uptake_time_delay <- merge(uptake_time_delay,uptake_time_rna,by='date_num',all.x = T,all.y=F)
  names(uptake_time_delay)
  
  # remove additional columns
  uptake_time_delay$dummy    <- NULL
  uptake_time_delay$date_num <- NULL
  
  # set imputed NA's as 0
  uptake_time_delay[is.na(uptake_time_delay)] <- 0
  
  return(uptake_time_delay)
}

get_vaccine_protection_matrix <- function(vaccine_uptake, dose_tag, parms, ndays_sim = NULL){
  
  if(dose_tag == '_A_'){
    delay_protection_rna   <- parms['delay_protection_rna1']    
    delay_protection_adeno <- parms['delay_protection_adeno1']
  } else if(dose_tag == '_B_'){
    delay_protection_rna   <- parms['delay_protection_rna2']    
    delay_protection_adeno <- parms['delay_protection_adeno2']
  } else if(dose_tag == '_E_'){
    delay_protection_rna   <- parms['delay_protection_rna2']    
    delay_protection_adeno <- parms['delay_protection_adeno2']
  } else if(dose_tag == '_F_'){
    delay_protection_rna   <- parms['delay_protection_rna2']    
    delay_protection_adeno <- parms['delay_protection_adeno2']
  } else{
    stop("VACCINE UPTAKE DOSE UNKNOWN:",dose_tag)
  }

  # check of given dose_tag is present, else, use dummy and set all info to 0
  if(!any(grepl(dose_tag,names(vaccine_uptake)))){
    dose_tag <- '_A_'
    vaccine_uptake[,grepl(dose_tag,names(vaccine_uptake))] <- 0
  }
  
  # protection (after uptake)
  vaccine_protection <- uptake2protection(vaccine_uptake[,grepl('date',names(vaccine_uptake)) | grepl(dose_tag,names(vaccine_uptake))],
                                          delay_protection_rna,
                                          delay_protection_adeno)
  # add time steps ahead of the given vaccine uptake.
  V_mat <- rbind(matrix(0, nrow = as.numeric(vaccine_protection$date[1]-get_scm_start_date()-1)/parms['h'], ncol = ncol(vaccine_protection[,-1])),
                 as.matrix(vaccine_protection[,-1]))
  
  
  # add time steps after the given vaccine uptake
  if(!is.null(ndays_sim) && (ndays_sim/parms['h']) > nrow(V_mat)){
    V_mat <- rbind(V_mat,
                   matrix(0, nrow = ceiling((ndays_sim/parms['h']) - nrow(V_mat)), ncol = ncol(V_mat)))
  }
  
  
  return(V_mat)
}

plot_uptake_by_age <- function(scen_out_f,
                               endpoint_reported_data = NA,
                               legend_label = 'Uptake'){
  
  age_opt     <- get_age_groups()
  ref_age_col <- paste0('vaccine_A_ages',age_opt)
  y_lim       <- c(0,1.5e6)
  x_ticks     <- pretty(range(scen_out_f$date),10)
  
  perc_breaks    <- pretty(0:100,7)
  y_ticks_perc   <- paste0(perc_breaks,'%')
  par(mar=c(5,5,4,5),mfrow=c(2,2))
  
  ref_data_be <- get_observed_incidence_data()
  
  i_age <- 5
  for(i_age in 10:3){
    
    y_ticks      <- perc_breaks * pop_data_be[i_age] / 100
    y_ticks_k    <- paste0(round(y_ticks/1e3),'k')
    y_ticks_k[1] <- 0
    y_ticks_k
    
    plot(ref_data_be$date+1,
         t(cumsum(ref_data_be[ref_age_col[i_age]])),
         ylab='Uptake',
         xlab='',
         xlim= range(scen_out_f$date),
         ylim = range(y_ticks),
         xaxt='n',
         yaxt='n',
         main=paste0(gsub('_','-',age_opt[i_age]),'y'))
    lines(scen_out_f$date+1,
          cumsum(scen_out_f[,11+i_age]),
          col=4,
          lwd=2)
    axis(1,x_ticks,format(x_ticks,'%d/%m'))
    axis(2,y_ticks,y_ticks_perc,las=2)
    axis(4,y_ticks,y_ticks_k,las=2)
    abline(v=x_ticks,lty=3,col='grey')
    abline(h=y_ticks,lty=3,col='grey')
    
    legend('topleft',
           c('Reported uptake',
             paste('SCM:', tolower(legend_label))),
           pch=c(1,NA),
           lwd=c(NA,2),
           col=c(1,4),
           cex=0.6,
           bg='white'
    )
    
    if(!is.na(endpoint_reported_data)){
      abline(v=endpoint_reported_data,lwd=2,col=3)  
      text(x = endpoint_reported_data,
           y = max(y_ticks)*0.98,
           'reported',
           pos=2,
           col=3,
           cex=0.7)
      text(x = endpoint_reported_data,
           y = max(y_ticks)*0.98,
           'scenario',
           pos=4,
           col=3,
           cex=0.7)
    }
  }
}

## VSC functions ----
is_vsc <- function(){
  return(nchar(system('echo $VSC_HOME',intern = T))>0)
}

get_vsc_scratch_folder <- function(){
  return(system('echo $VSC_SCRATCH',intern = T))
}

is_vsc_vaughan_cluster <- function(){
  return(system('echo $VSC_INSTITUTE_CLUSTER',intern = T) == "vaughan")
}

get_output_dir <- function(...){
  
  # define output directory name
  if(is_vsc_vaughan_cluster()){
    output_dir <- file.path(get_vsc_scratch_folder(),'stochastic_model/output',...,'')
  } else{
    output_dir <- file.path('./output',...,'')
  }
  
  # return the output directory name
  return(output_dir)
}

init_output_dir <- function(...){
  
  # get output directory name
  output_dir <- get_output_dir(...)
  
  # create directory if not existing yet (recursive)
  if(!dir.exists(output_dir)){
    dir.create(output_dir,recursive = T)
  }
  
  # return the output directory name
  return(output_dir)
}

## WGS functions ----
approx_sequenced <- function(date,sequenced_week){
  
  out <- approx(x = date[!is.na(sequenced_week)],
                y = sequenced_week[!is.na(sequenced_week)],
                xout = seq(min(date),max(date),1))
  
  return(data.frame(x = as.character(out$x),
                    y = round(out$y)))
}

## HOSPITAL ADMISSION functions ----
# hospital admission rate: DEFAULT
# ==>> (1 - exp(-h*(1-phi0)*delta2)
#
# hospital admission rate: WITH VOC:
# ==>> (1 - exp(-h*((1-phi0) + phi0*VOC_hosp)*delta2)
#
# Odds ratio for hospital admission: ONE VOC
#     (1 - exp(-h*((1-phi0) + phi0*VOC_new)*delta2)) /
#     (1 - exp(-h*(1-phi0)*delta2))
#
# Odds ratio for hospital admission: TWO VOCs
#     (1 - exp(-h*((1-phi0) + phi0*VOC_new)*delta2)) /
#     (1 - exp(-h*((1-phi0) + phi0*VOC_previous)*delta2))
get_log_VOC_hosp <- function(phi0,delta2,h,odds_ratio,log_VOC_previous=NA){
  
  # defensive check
  if(length(odds_ratio) > 1 && length(odds_ratio) != length(phi0)){
    stop("THE PROVIDED HOSPITAL ADMISSION ODDS RATIO IS NOT VALID")
  }
  
  # define optional secondary VOC factor
  VOC_previous <- ifelse(is.na(log_VOC_previous),0,expit(log_VOC_previous))
  
  # get VOC hospital admission factor
  VOC_new <- log(1- odds_ratio * (1 - exp(-h*((1-phi0)+phi0*VOC_previous)*delta2)))/(-h*delta2*phi0) - (1-phi0)/phi0
  
  # use cutoff of 1e-15, to prevent NaN with logit
  VOC_new[VOC_new <= 1e-15] <- 1e-15
  
  # convert to logit
  log_VOC_new <- logit(VOC_new)
  
  return(log_VOC_new)
}

# help function
get_VOC_hosp <- function(phi0,delta2,h,odds_ratio,log_VOC_previous=NA){
  return(expit(get_log_VOC_hosp(phi0,delta2,h,odds_ratio,log_VOC_previous)))
}
  
# see comments with "get_log_VOC_hosp"
get_VOC_hosp_odds_ratio <- function(phi0,delta2,h,log_VOC_new,log_VOC_previous=NA){
  
  VOC_new       <- expit(log_VOC_new)
  
  if(any(is.na(log_VOC_previous))){
    VOC_previous <- 0
  } else {
    VOC_previous <- expit(log_VOC_previous)
  }
  
  OR <- as.numeric((1 - exp(-h*((1-phi0) + phi0*VOC_new)*delta2)) /
                     (1 - exp(-h*((1-phi0) + phi0*VOC_previous)*delta2)))
  
  return(OR)
}

## LOG FUNCTIONS ----
write_log_file <- function(log_file_name, file_append = FALSE) {
    
  # load package info
  pkg_loaded         <- (.packages())
  pkg_all_version    <- installed.packages()[,'Version']
  pkg_loaded_version <- pkg_all_version[names(pkg_all_version) %in% pkg_loaded]
  
  # load git info
  
  
  # combine into strings
  pkg_txt               <- paste0(names(pkg_loaded_version),':\t', pkg_loaded_version)
  pkg_nchar             <- unlist(lapply(names(pkg_loaded_version),nchar))
  pkg_txt[pkg_nchar<=5] <- gsub(':\t',':\t\t',pkg_txt[pkg_nchar<=5])
  
  # write to file
  cat(R.Version()$version.string,file=log_file_name,append = file_append, fill = T)
  cat(paste0('platform:\t',R.Version()$platform),file=log_file_name,append = TRUE,fill = T)
  cat(paste0('current date:\t',Sys.Date()),file=log_file_name,append = TRUE,fill = T)
  cat('\nLOADED PACKAGES',file=log_file_name,append = TRUE,fill = T)
  sapply(pkg_txt,FUN=cat,file=log_file_name,append = TRUE,fill=TRUE) -> tmp
}


# COMPARTMENT AND TRANSITION FUNCTIONS ----
get_basic_compartment_names <- function(col_tag="",bool_incidence=TRUE){
  
  # default
  icol_names_ode <- c('S', 'E', 'I_presym', 'I_asym', 'I_mild', 
                  'I_sev', 'I_hosp', 'I_icu', 'D', 'R')
  
  # Additional compartments to track incidence
  icol_names_incidence <- c('new_hosp_total', 'asym_mild_infected', 
                            'hosp_recovered')
  
  if(nchar(col_tag)>0){
    #TODO: make uniform
    icol_names_ode <- paste0(substr(icol_names_ode,0,1),col_tag,substr(icol_names_ode,2,100))
    
    # account for col_tag
    icol_names_incidence <- paste(icol_names_incidence,col_tag,sep='_')
  }
  
  # join
  if(bool_incidence){
    icol_names <- c(icol_names_ode,icol_names_incidence)
  } else{
    icol_names <- c(icol_names_ode)
  }
  
  # return
  return(icol_names)
}

get_all_compartment_names <- function(num_age_groups = 10,bool_incidence = TRUE){
  icol_names_all <- c(get_basic_compartment_names(bool_incidence=TRUE),
                      get_basic_compartment_names(col_tag='voc',bool_incidence=TRUE)[-1])     # remove Svoc
  
  if(bool_incidence){
    icol_names_all <- c(icol_names_all,
                        paste0('mild_infected',c('','_voc')),
                        paste0('recov_mild_infection',c('','_voc')),
                        paste0('new_icu',c('','_voc')),
                        paste0('new_D',c('','_voc')),
                        paste0('new_E',c('','_voc')))
  }
  
  return(icol_names_all)
}

get_transition_names <- function(){
  
  icol_names <- get_all_compartment_names()
  
  t_names <- c(icol_names[2:8],
               paste(icol_names[9],c('hosp','icu'),sep='_'),
               paste(icol_names[10],c('asym','mild','hosp','icu'),sep='_'))
  t_names <- paste0('t_',t_names)
  t_names <- c(t_names,paste0(substr(t_names,0,3),'voc',substr(t_names,4,20)))
  
  return(t_names)
}

get_compartment_index <- function(num_age_groups = 10){
  icol_names_all <- get_all_compartment_names()

  icol_tmp_vec <- matrix(1:(length(icol_names_all)*num_age_groups),
                         ncol=length(icol_names_all),
                         nrow=num_age_groups,
                         byrow = F)
  
  icol_tmp_vec<- data.frame(icol_tmp_vec)
  names(icol_tmp_vec) <- icol_names_all
  return(icol_tmp_vec)
}

# to distinguish the core ODE compartments from the ones to track incidence
is_core_ODE <- function(icol_tmp_vec){
  
  # if first character is a capital letter, it is part of the basic ODE
  # else, it is a compartment to track incidence
  first_char <- substr(names(icol_tmp_vec),0,1)
  bool_ode  <- first_char == toupper(first_char)
  
  return(bool_ode)
  
  }

get_population_matrix <- function(num_age_groups=10){
  compartment_names  <- get_all_compartment_names()
  pop_matrix         <- matrix(0,nrow=num_age_groups,ncol=length(compartment_names))
  colnames(pop_matrix)  <- compartment_names

  return(pop_matrix)
}

get_transition_matrix <- function(num_age_groups=10){
  transition_names    <- get_transition_names()
  transition_matrix   <- data.frame(matrix(0,nrow=num_age_groups,ncol=length(transition_names)))
  names(transition_matrix) <- transition_names
  
  return(transition_matrix)
}

update_compartments <- function(pop_state,new_transitions,num_age_groups = 10){
  
  pop_matrix         <- get_population_matrix(pop_state,num_age_groups=num_age_groups)
  transition_matrix  <- get_transition_matrix(new_transitions,num_age_groups=num_age_groups)

  pop_matrix <- update_compartments_matrix(pop_matrix,transition_matrix)
  
  # adjust to the original structure
  pop_state <- unlist(pop_matrix)
  
  # return
  return(pop_state)
}

# note: we make use of matrix types (<double>) instead of data.frame (<lists>)
# to speed up the code. The disadvantage is that we cannot use the '$' operator
# but have to explicitly call the name of the compartment, which is error prone.
#transition_matrix <- transition_array[,,i_pop_type]
update_compartments_matrix <- function(pop_matrix,transition_matrix){
    
  # update the population, by starting with empty matrix
  pop_matrix_new          = pop_matrix*0
  
  # account for non-VOC transitions
  pop_matrix_new[,'S']        = - transition_matrix[,'t_E'] - transition_matrix[,'t_Evoc']
  pop_matrix_new[,'E']        = transition_matrix[,'t_E'] - transition_matrix[,'t_I_presym']
  pop_matrix_new[,'I_presym'] = transition_matrix[,'t_I_presym'] - transition_matrix[,'t_I_asym'] - transition_matrix[,'t_I_mild']
  pop_matrix_new[,'I_asym']   = transition_matrix[,'t_I_asym'] - transition_matrix[,'t_R_asym']
  pop_matrix_new[,'I_mild']   = transition_matrix[,'t_I_mild'] - transition_matrix[,'t_I_sev']  - transition_matrix[,'t_R_mild']
  pop_matrix_new[,'I_sev']    = transition_matrix[,'t_I_sev']  - transition_matrix[,'t_I_hosp'] - transition_matrix[,'t_I_icu']
  pop_matrix_new[,'I_hosp']   = transition_matrix[,'t_I_hosp'] - transition_matrix[,'t_R_hosp'] - transition_matrix[,'t_D_hosp']
  pop_matrix_new[,'I_icu']    = transition_matrix[,'t_I_icu']  - transition_matrix[,'t_R_icu']  - transition_matrix[,'t_D_icu']
  pop_matrix_new[,'R']        = transition_matrix[,'t_R_asym'] + transition_matrix[,'t_R_mild'] + transition_matrix[,'t_R_hosp'] #+ transition_matrix[,'t_R_icu # TODO !!
  pop_matrix_new[,'D']        = transition_matrix[,'t_D_hosp'] + transition_matrix[,'t_D_icu']  
  
  # account for VOC transitions
  pop_matrix_new[,'Evoc']        = transition_matrix[,'t_Evoc'] - transition_matrix[,'t_Ivoc_presym']
  pop_matrix_new[,'Ivoc_presym'] = transition_matrix[,'t_Ivoc_presym'] - transition_matrix[,'t_Ivoc_asym'] - transition_matrix[,'t_Ivoc_mild']
  pop_matrix_new[,'Ivoc_asym']   = transition_matrix[,'t_Ivoc_asym'] - transition_matrix[,'t_Rvoc_asym']
  pop_matrix_new[,'Ivoc_mild']   = transition_matrix[,'t_Ivoc_mild'] - transition_matrix[,'t_Ivoc_sev']  - transition_matrix[,'t_Rvoc_mild']
  pop_matrix_new[,'Ivoc_sev']    = transition_matrix[,'t_Ivoc_sev']  - transition_matrix[,'t_Ivoc_hosp'] - transition_matrix[,'t_Ivoc_icu']
  pop_matrix_new[,'Ivoc_hosp']   = transition_matrix[,'t_Ivoc_hosp'] - transition_matrix[,'t_Rvoc_hosp'] - transition_matrix[,'t_Dvoc_hosp']
  pop_matrix_new[,'Ivoc_icu']    = transition_matrix[,'t_Ivoc_icu']  - transition_matrix[,'t_Rvoc_icu']  - transition_matrix[,'t_Dvoc_icu']
  pop_matrix_new[,'Rvoc']        = transition_matrix[,'t_Rvoc_asym'] + transition_matrix[,'t_Rvoc_mild'] + transition_matrix[,'t_Rvoc_hosp'] #+ transition_matrix[,'t_Rvoc_icu # TODO !!
  pop_matrix_new[,'Dvoc']        = transition_matrix[,'t_Dvoc_hosp'] + transition_matrix[,'t_Dvoc_icu']  
  
  # update population matrix
  pop_matrix <- pop_matrix + pop_matrix_new
  
  # track incidence 
  pop_matrix[,'new_hosp_total']     <- transition_matrix[,'t_I_hosp'] + transition_matrix[,'t_I_icu']
  pop_matrix[,'asym_mild_infected'] <- transition_matrix[,'t_I_asym'] + transition_matrix[,'t_I_mild']
  pop_matrix[,'hosp_recovered']     <- transition_matrix[,'t_R_hosp'] + transition_matrix[,'t_R_icu']
  
  pop_matrix[,'new_hosp_total_voc']     <- transition_matrix[,'t_Ivoc_hosp'] + transition_matrix[,'t_Ivoc_icu']
  pop_matrix[,'asym_mild_infected_voc'] <- transition_matrix[,'t_Ivoc_asym'] + transition_matrix[,'t_Ivoc_mild']
  pop_matrix[,'hosp_recovered_voc']     <- transition_matrix[,'t_Rvoc_hosp'] + transition_matrix[,'t_Rvoc_icu']
  
  # track incidence 
  pop_matrix[,'mild_infected']            <- transition_matrix[,'t_I_mild']
  pop_matrix[,'mild_infected_voc']        <- transition_matrix[,'t_Ivoc_mild']
  pop_matrix[,'recov_mild_infection']     <- transition_matrix[,'t_R_mild']
  pop_matrix[,'recov_mild_infection_voc'] <- transition_matrix[,'t_Rvoc_mild']
  pop_matrix[,'new_icu']                  <- transition_matrix[,'t_I_icu']
  pop_matrix[,'new_icu_voc']              <- transition_matrix[,'t_Ivoc_icu']
  
  # track mortality
  pop_matrix[,'new_D']                    <- pop_matrix_new[,'D'] 
  pop_matrix[,'new_D_voc']                <- pop_matrix_new[,'Dvoc'] 
  
  # track infections
  pop_matrix[,'new_E']                    <- transition_matrix[,'t_E']
  pop_matrix[,'new_E_voc']                <- transition_matrix[,'t_Evoc']

  # return
  return(pop_matrix)
}


get_3dim_nfull_array <- function(nfull_array_4dim){
  if(length(dim(nfull_array_4dim))==4){
    nfull_array_3dim <- nfull_array_4dim[,,,1]
    if(dim(nfull_array_4dim)[4]>1){
      for(i_vac in 2:dim(nfull_array_4dim)[4]){
        nfull_array_3dim <- nfull_array_3dim + nfull_array_4dim[,,,i_vac]
      }  
    }
    return(nfull_array_3dim)
  } else{
    return(nfull_array_4dim) # already 3 dimensions
  }
}

get_cases_summary <- function(nfull_array_aggr,state_name,bool_extend = TRUE,col_prefix=''){

  # if multiple vaccine-related health classes are present, aggregate
  nfull_array_aggr <- get_3dim_nfull_array(nfull_array_aggr)

  if(bool_extend){
    state_name <- c(state_name,paste0(state_name,'_voc'))
    state_name <- c(state_name,paste0(substr(state_name,0,1),'voc',substr(state_name,2,100)))
    state_name <- state_name[state_name %in% get_all_compartment_names()]
  }

  matrix_time_age <- nfull_array_aggr[,,state_name[1]]
  if(length(state_name)>1){
    for(i_state in state_name[-1]){
      matrix_time_age <- matrix_time_age + nfull_array_aggr[,,i_state]
    }
   }

  nsim_day_id <- 1:nrow(nfull_array_aggr)-1
  df_time_age <- data.frame(day=nsim_day_id,cases=matrix_time_age)
  
  if(nchar(col_prefix)>0){
    names(df_time_age) <- gsub('cases',col_prefix,names(df_time_age))
  }

  return(df_time_age)
}


get_transition_origin <- function(){
  ## Transition probabilities
  col_vec_origin = c('S',
                'E',

                # Ipresym to I_asym en I_mild
                'I_presym',
                'I_presym',

                # Severe infections
                'I_mild',

                # Hospital and ICU admission
                'I_sev', 
                'I_sev',

                # Mortality
                'I_hosp',
                'I_icu', 

                # Recovery
                'I_asym',
                'I_mild',
                'I_hosp',
                'I_icu',

                # idem for VoC
                'S',                   
                'Evoc','Ivoc_presym',   
                'Ivoc_presym',          
                'Ivoc_mild','Ivoc_sev', 
                'Ivoc_sev','Ivoc_hosp','Ivoc_icu',
                'Ivoc_asym', 'Ivoc_mild', 
                'Ivoc_hosp', 'Ivoc_icu'
  );
  return(col_vec_origin)
}

# Hospital probability 
# phi0 # age-specific probability NOT to be hospitalized
get_phi0_vac <- function(phi0,ve_hosp){
  return(phi0 + (1-phi0)*ve_hosp)
}

#Sensitivity function (serial serological data)
serology_sens_function <- function(t){
  return(expit(-8.397 + 0.928*t))
}  

get_prob_array <- function(parms,log_VOC_XXX_hosp,calendar_time,pop_categories,ve_tag='ve'){
  
  # subset parameter names
  parms_names <- names(parms)
  
  ## General model parameters 
  gamma  = exp(parms['log_gamma']);   # average length of the latency period - 2 days
  theta  = exp(parms['log_theta']);   # infectious period at pre-symptomatic stage - 3.5 days
  delta1 = exp(parms['log_delta1']);  # infectious period at asymptomatic stage - 3.5 days
  delta2 = exp(parms['log_delta2']);  # infectious period at (mild) symptomatic stage - 3.5 days
  
  phi0  = expit(parms[grepl('log_phi0',parms_names)]);  # proportion of symptomatically infected with mild symptoms
  p_vec = parms[grepl('p_asympt_age',parms_names)];     # proportion of asymptomatic cases
  h     = parms['h'];                                   # resolution of the binomial chains
  
  omega = exp(parms[grepl('log_omega',parms_names)]);          # waiting time between symptom onset and hospitalization
  phi1  = expit(parms[grepl('log_phi1_age',parms_names)] * parms['log_phi1_add']);  # proportion of severly infected with regular hospitalization
  
  mu_param = parms[grepl('mu_sev',parms_names)]
  mu_sev = c(0,expit(mu_param[c(1,1:8)]));   # age-specific hospital fatality ratio - probability of dying
  delta3 = exp(parms['log_delta3']);         # recovery period in severely infected/length of stay - 7 days
  
  sel_parms  = parms
  VOC_start  = parms[grepl('VOC_.*_start',parms_names)]
  
  # option to get VOC parameters
  if(calendar_time %in% VOC_start){
    VOC_name             <- gsub('_start','',names(VOC_start))
    sel_VOC              <- VOC_name[VOC_start == calendar_time]
    sel_VOC_param        <- parms[grepl(sel_VOC,parms_names)]
    names(sel_VOC_param) <- gsub(paste0(sel_VOC,'_'),'',names(sel_VOC_param))
    names(sel_VOC_param) <- gsub(paste0('_',sel_VOC),'',names(sel_VOC_param))
    sel_VOC_param
    
    # VOC parameters  
    VOC_hosp    <- expit(log_VOC_XXX_hosp[which(VOC_start == calendar_time),])
    phi0        <- phi0 - (phi0 * VOC_hosp)
    phi1        <- expit(parms[grepl('log_phi1_age',parms_names)] * sel_VOC_param['log_phi1_add']);  # proportion of severly infected with regular hospitalization
    delta3      <- exp(sel_VOC_param['log_delta3']);    # recovery rate from hospital (from delta VOC)
    gamma       <- exp(parms['log_gamma']) * ifelse('log_gamma_factor' %in% names(sel_VOC_param), exp(sel_VOC_param['log_gamma_factor']), 1)
    
    log_VOC_mu_param <- sel_VOC_param[grepl('log_mu_sev',names(sel_VOC_param))]
    mu_sev           <- c(0,expit(log_VOC_mu_param[c(1,1:8)])); # age-specific hospital fatality ratio - probability of dying
    
    sel_parms = sel_VOC_param
  } 
  
  # default transition probabilities/rates
  prob_matrix <- cbind(
    t_E        = NA,        
    t_I_presym =  (1 - exp(-h*gamma)),
    
    t_I_asym   =  (1 - exp(-h*p_vec*theta)),
    t_I_mild   =  (1 - exp(-h*(1-p_vec)*theta)),
    
    t_I_sev    =  (1 - exp(-h*(1-phi0)*delta2)),
    
    t_I_hosp   =  (1 - exp(-h*phi1*omega)),
    t_I_icu    =  (1 - exp(-h*(1 - phi1)*omega)),
    
    t_D_hosp   =  (1 - exp(-h*delta3*mu_sev)),
    t_D_icu    =  (1 - exp(-h*delta3*mu_sev)),
    
    t_R_asym   =  (1 - exp(-h*delta1)),
    t_R_mild   =  (1 - exp(-h*phi0*delta2)),
    t_R_hosp   =  (1 - exp(-h*delta3*(1-mu_sev))),
    t_R_icu    =  (1 - exp(-h*delta3*(1-mu_sev)))
  )
  
  # duplicate for vaccination states
  prob_array <- array(prob_matrix,dim=c(nrow(prob_matrix),ncol(prob_matrix),10),
                      dimnames = list(NULL,
                                      colnames(prob_matrix),
                                      pop_categories$name))
  
  # adjust for infection- or vaccine-related protection
  for(i_pop_catetory in pop_categories$name[-1]){
    param_name <- paste0(ve_tag,'_incr_severe_',gsub('vac_','',i_pop_catetory))
    prob_array[,'t_I_sev',i_pop_catetory]   <- 1 - exp(-h*(1-get_phi0_vac(phi0,sel_parms[param_name]))*delta2)
    prob_array[,'t_R_mild',i_pop_catetory]  <- 1 - exp(-h*   get_phi0_vac(phi0,sel_parms[param_name]) *delta2)
  }
  
  # return
  return(prob_array)
}

sum_by_pop_type <- function(pop_array){
  return(apply(pop_array,2,rowSums))
}



