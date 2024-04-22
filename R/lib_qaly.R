########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is to CALCULATE QALY output
#
# Copyright 2024, SIMID                                       
########################################################################### #

if('EQ5D.be' %in% installed.packages()){
  library(EQ5D.be) 
} else{
  source('data/eq5d/eq5d_popnormINDEX.R') #backup script for older platforms
}

multi_qaly_loss <- function(output_files, 
                            col_lib,
                            i_row_reference=2,
                            time_horizon,
                            file_tag,
                            axis_lim = NA,
                            bool_plot = TRUE){
  

  # reference QALY loss
  reference_qaly_loss       <- get_qaly_loss(output_files =  output_files,
                                             scen_tag = col_lib$tags[i_row_reference],
                                             time_horizon = time_horizon)
  # initialise aggregated table
  aggregated_qaly_gain <- NULL
  
  i_col <- 2
  for(i_col in 2:nrow(col_lib)){
      scenario_qaly_loss <- get_qaly_loss(output_files =  output_files,
                                          scen_tag = col_lib$tags[i_col],
                                          time_horizon=time_horizon)
      
      # calculate difference
      scenario_qaly_gain <- reference_qaly_loss - scenario_qaly_loss
      
      # add labels et
      scenario_qaly_gain$label <- col_lib$label[i_col]
      
      # row bind
      aggregated_qaly_gain    <- rbind(aggregated_qaly_gain,scenario_qaly_gain)
  }
  
  scen_order <- as.numeric(factor(unique(aggregated_qaly_gain$label),levels= col_lib$label[col_lib$label %in% unique(aggregated_qaly_gain$label)],ordered=TRUE))
  if(any(is.na(axis_lim))){
    axis_lim <- range(aggregated_qaly_gain$total)
  }
  
  # summary box-plot figures
  if(bool_plot){
    par(mar=c(5,19,1,1),mfrow=c(1,1))  
    bplot <- boxplot(aggregated_qaly_gain$total ~ aggregated_qaly_gain$label,
            las=1,horizontal = T,xlab='QALY gain',ylab='',
            cex.axis=0.8,ylim=axis_lim)  
    grid(nx=NULL,ny=NA)
    abline(v=0,lty=2)
    aggr_mean <- aggregate(total ~ label,data=aggregated_qaly_gain,mean)
    if(all(bplot$names == aggr_mean$label)){
      points(x=aggr_mean$total,y=nrow(aggr_mean):1,pch=4)
    }
    
    boxplot(aggregated_qaly_gain$total_disc ~ aggregated_qaly_gain$label,
            las=1,horizontal = T,xlab='QALY gain (discounted)',ylab='',
            at=rev(scen_order),cex.axis=0.8)  
    grid(nx=NULL,ny=NA)
    abline(v=0,lty=2)
    aggr_mean <- aggregate(total ~ label,data=aggregated_qaly_gain,mean)
    if(all(bplot$names == aggr_mean$label)){
      points(x=aggr_mean$total,y=nrow(aggr_mean):1,pch=4)
    }
  }
  
  # add burden of disease
  reference_qaly_loss$label <- 'Estimated BoD'
  aggregated_qaly_gain <- rbind(aggregated_qaly_gain,reference_qaly_loss)
  scen_order <- c(scen_order,length(scen_order)+1)
  
  # plot (prevented) burden per scenario
  if(bool_plot){
    
    # burden by age
    par(mar=c(5,5,1,1),mfrow=c(3,1)) 
    
    sel_columns_age  <- names(aggregated_qaly_gain)[grepl('qaly.*age',names(aggregated_qaly_gain)) & !grepl('all',names(aggregated_qaly_gain))]
    sel_columns_qaly_illness <- sel_columns_age[!grepl('mort',sel_columns_age)]
    sel_columns_qaly_mort    <- sel_columns_age[grepl('mort',sel_columns_age)]
    sel_columns_num_mort     <- names(aggregated_qaly_gain)[grepl('n_mort_age',names(aggregated_qaly_gain))]
    
    names_matrix <- matrix(sel_columns_qaly_illness,ncol=10,byrow = T)
    bar_labels <- gsub('_age.*','',names_matrix[,1])
    bar_labels <- gsub('qaly_*','',bar_labels)
    
    aggregated_qaly_gain_mean  <- aggregate(. ~ label, data = aggregated_qaly_gain,mean)
    aggregated_qaly_gain_CI_LL <- aggregate(. ~ label, data = aggregated_qaly_gain,quantile,0.025)
    aggregated_qaly_gain_CI_UL <- aggregate(. ~ label, data = aggregated_qaly_gain,quantile,0.975)
    
    opt_label <- aggregated_qaly_gain_mean$label
    
    get_qaly_gain_matrix_summary <- function(aggregated_qaly_gain_mean){
      
      qaly_gain_matrix_summary <- NULL
      
      i_row <- 1
      for(i_row in 1:nrow(aggregated_qaly_gain_mean)){
        qaly_gain_matrix <- matrix(unlist(aggregated_qaly_gain_mean[i_row,sel_columns_qaly_illness]),ncol=10,byrow = T)
        colnames(qaly_gain_matrix) <- get_age_groups()
        qaly_gain_sum <- colSums(qaly_gain_matrix)
        
        num_mort_prevented_matrix <- unlist(aggregated_qaly_gain_mean[i_row,sel_columns_num_mort])
        names(num_mort_prevented_matrix) <- get_age_groups()
        
        qaly_mort_prevented_matrix <- unlist(aggregated_qaly_gain_mean[i_row,sel_columns_qaly_mort])
        names(qaly_mort_prevented_matrix) <- get_age_groups()
        
        qaly_gain_matrix_summary<- rbind(qaly_gain_matrix_summary,
                                         data.frame(scenario = opt_label[i_row],
                                                    burden = c(bar_labels,'total','num_deaths','qaly_mort'),
                                                    age = rbind(qaly_gain_matrix,
                                                                qaly_gain_sum,
                                                                num_mort_prevented_matrix,
                                                                qaly_mort_prevented_matrix)))
        
      }
      return(qaly_gain_matrix_summary)
    }
    
    qaly_gain_matrix_summary <- get_qaly_gain_matrix_summary(aggregated_qaly_gain_mean)
    qaly_gain_matrix_summary_CI_LL <- get_qaly_gain_matrix_summary(aggregated_qaly_gain_CI_LL)
    qaly_gain_matrix_summary_CI_UL <- get_qaly_gain_matrix_summary(aggregated_qaly_gain_CI_UL)
    
    y_lim_qaly <- range(qaly_gain_matrix_summary[qaly_gain_matrix_summary$scenario != 'Estimated BoD' & qaly_gain_matrix_summary$burden == 'total',-(1:2)])
    y_lim_mort <- range(qaly_gain_matrix_summary[qaly_gain_matrix_summary$scenario != 'Estimated BoD' & qaly_gain_matrix_summary$burden == 'num_deaths',-(1:2)])
    y_lim_qaly_mort <- range(qaly_gain_matrix_summary[qaly_gain_matrix_summary$scenario != 'Estimated BoD' & qaly_gain_matrix_summary$burden == 'qaly_mort',-(1:2)])
    
    i_row = 1
    for(i_row in 1:length(opt_label)){
      
      qaly_gain_matrix <- qaly_gain_matrix_summary[qaly_gain_matrix_summary$scenario == opt_label[i_row] &
                                                     qaly_gain_matrix_summary$burden %in% bar_labels,-(1:2)]
      
      qaly_gain_matrix <- matrix(unlist(qaly_gain_matrix),ncol=10,byrow = F)
      colnames(qaly_gain_matrix) <- get_age_groups()
      
      if(grepl('Increased',opt_label[i_row])){
        y_lim <- c(0,max(y_lim_qaly))
      } else{
        y_lim <- NULL
      }
      
      # stacked
      barplot(qaly_gain_matrix,
              xlab = 'age group',
              ylab = ifelse(i_row == 1,'QALY loss','QALY gain'),
              main=opt_label[i_row],
              cex.axis=0.8,
              ylim = y_lim)
      grid(nx=NA,ny=NULL)
      legend('topleft',rev(bar_labels),fill=rev(grey.colors(length(bar_labels))))
      
      # quality adjusted life years lost
      sel_summary <- qaly_gain_matrix_summary$scenario == opt_label[i_row] & qaly_gain_matrix_summary$burden %in% 'qaly_mort'
      sel_label <- opt_label[i_row]
      y_lab <- 'debug'
      
      add_barplot_whiskers <- function(barplot_ticks,bar_values_ul,bar_values_ll){
        
        bar_values_ll <- unlist(bar_values_ll)
        bar_values_ul <- unlist(bar_values_ul)
        
        is_valid <- (bar_values_ll - bar_values_ul) != 0
        if(any(is_valid)){
          arrows(barplot_ticks[is_valid],unlist(bar_values_ul[is_valid]), 
                 barplot_ticks[is_valid],unlist(bar_values_ll[is_valid]), 
                 lwd = 1, angle = 90,
                 code = 3, length = 0.05)          
        }
      }
      
      plot_bar_CI <- function(sel_summary,sel_label,y_lab){
        
        bar_values <- unlist(qaly_gain_matrix_summary[sel_summary,-(1:2)])
        
        # error-bars
        bar_values_ll <- unlist(qaly_gain_matrix_summary_CI_LL[sel_summary,-(1:2)])
        bar_values_ul <- unlist(qaly_gain_matrix_summary_CI_UL[sel_summary,-(1:2)])
        
        # add group names
        names(bar_values) <- get_age_groups()
        
        y_lim <- range(c(0,bar_values,bar_values_ll,bar_values_ul))
        
        # plot bars
        bplt <- barplot(bar_values,
                        xlab = 'age group',
                        # ylab = ifelse(i_row == 1,'Mortality-related QALY lost','Mortality-related QALY gained'),
                        ylab = y_lab,
                        main = sel_label[i_row],
                        ylim = y_lim)
        grid(nx=NA,ny=NULL)
        
        add_barplot_whiskers(barplot_ticks = bplt,
                             bar_values_ll = bar_values_ll,
                             bar_values_ul = bar_values_ul)
      }
      
      # number of deaths
      sel_summary <- qaly_gain_matrix_summary$scenario == opt_label[i_row] & qaly_gain_matrix_summary$burden %in% 'num_deaths'
      plot_bar_CI(sel_summary,sel_label = opt_label[i_row],y_lab = ifelse(i_row == 1,'Number of deaths','Prevented deaths'))
      
      
      # quality adjusted life years lost
      sel_summary <- qaly_gain_matrix_summary$scenario == opt_label[i_row] & qaly_gain_matrix_summary$burden %in% 'qaly_mort'
      plot_bar_CI(sel_summary,sel_label = opt_label[i_row],y_lab = ifelse(i_row == 1,'Mortality-related QALY lost','Mortality-related QALY gained'))
      
    }
    
    qaly_gain_df <- qaly_gain_matrix_summary[grepl('increased',tolower(qaly_gain_matrix_summary$scenario)) &
                                                            qaly_gain_matrix_summary$burden=='total',]
    
    if(nrow(qaly_gain_df)>0){
      par(mfrow=c(3,1))
      
      qaly_gain_matrix <- t(t(qaly_gain_df[,-(1:2)]))
      colnames(qaly_gain_matrix) <- get_age_group_labels()
      rownames(qaly_gain_matrix) <- qaly_gain_df$scenario
      barplot(qaly_gain_matrix,
              beside = T,
              xlab = 'Age group',
              ylab = 'QALY gain')
      grid(nx=NA,ny=NULL)
      legend('topright',qaly_gain_df$scenario,fill=(grey.colors(nrow(qaly_gain_df))),cex=0.7)
      
      # number of deaths averted
      num_mort_prevented_df <- qaly_gain_matrix_summary[grepl('increased',tolower(qaly_gain_matrix_summary$scenario)) &
                                                      qaly_gain_matrix_summary$burden=='num_deaths',]
      mort_prevented_matrix <- t(t(num_mort_prevented_df[,-(1:2)]))
      colnames(mort_prevented_matrix) <- get_age_groups()
      barplot(mort_prevented_matrix,
              beside = T,
              xlab = 'Age group',
              ylab = ifelse(i_row == 1,'Number of deaths','Prevented deaths'))
      grid(nx=NA,ny=NULL)
      legend('topleft',qaly_gain_df$scenario,fill=(grey.colors(nrow(qaly_gain_df))),cex=0.7)
    
      # quality adjusted life years lost
      sel_summary <- grepl('increased',tolower(qaly_gain_matrix_summary$scenario)) & qaly_gain_matrix_summary$burden=='qaly_mort'
      mort_qaly_gain_df <- qaly_gain_matrix_summary[sel_summary,-(1:2)]
      mort_qaly_gain_df <- t(t(mort_qaly_gain_df))
      colnames(mort_qaly_gain_df) <- get_age_group_labels()
  
      y_lim <- c(min(mort_qaly_gain_df),max(pretty(mort_qaly_gain_df)))
      
      bplot <- barplot(mort_qaly_gain_df,
                      beside = T,
                      xlab = 'Age group',
                      ylab = ifelse(i_row == 1,'Mortality-related QALY lost','Mortality-related QALY gain'),
                      ylim = y_lim)
      grid(nx=NA,ny=NULL)
      
      # incremental effect (illness and mortality) ----
      par(mfrow=c(1,1))
      qaly_out <- matrix(0,ncol=29,nrow=4)
      
      i_scen                 <- seq(1,ncol(qaly_gain_matrix)*3,3)
      qaly_out[1:2,i_scen]   <- rbind(mort_qaly_gain_df[1,],qaly_gain_matrix[1,])
      qaly_out[3:4,i_scen+1] <- rbind(mort_qaly_gain_df[2,],qaly_gain_matrix[2,])

      y_lim <- c(0,max(pretty(colSums(qaly_out),na.rm=T)))
      
      bplot_qaly <- barplot(qaly_out,beside = F,
                            col=c('darkred','darkred','lightblue','lightblue','white','white'),
                            density = c(200,30,200,30,0,0),
                            #border = c('darkred','darkred','lightblue','lightblue',NA,NA),
                            #border = NA,
                            xlab='Age group',
                            ylab='QALY gain',
                            ylim = y_lim,
                            las=2,
                            xpd=FALSE)
      grid(nx = NA,ny=NULL)
      abline(h=0,lwd=2)
      axis(1,at=bplot_qaly[i_scen]+0.5,labels=get_age_group_labels(),las=1,cex.axis=1.0)

      legend('topleft',
             paste0(rep(qaly_gain_df$scenario,each=2),
               c(': illness-related QALY',
               ': mortality-related QALY')),
             fill=c('darkred','darkred','lightblue','lightblue'),
             density=c(30,200,30,200),
             cex=0.8)
      
      # TOTAL
      bplot_qaly <- barplot(colSums(qaly_out),beside = F,
                            col=c('darkred','lightblue','white'),
                            #density = c(100,30,100,30,0,0),
                            #border = c('darkred','darkred','lightblue','lightblue',NA,NA),
                            #border = NA,
                            xlab='Age group',
                            ylab='QALY gain',
                            ylim = y_lim,
                            las=2,
                            xpd=FALSE)
      
      grid(nx = NA,ny=NULL)
      abline(h=0,lwd=2)
      axis(1,at=bplot_qaly[i_scen]+0.5,labels=get_age_group_labels(),las=1,cex.axis=1.0)
      
      legend('topleft',
             qaly_gain_df$scenario,
             fill=c('darkred','lightblue'),
             cex=0.8)
      
    }
  }
  
  # numeric results
  qaly_gain_mean    <- aggregate(. ~ label, data = aggregated_qaly_gain, mean)[scen_order,]
  qaly_gain_CrI_min <- aggregate(. ~ label, data = aggregated_qaly_gain, quantile,0.025)[scen_order,]
  qaly_gain_CrI_max <- aggregate(. ~ label, data = aggregated_qaly_gain, quantile,0.975)[scen_order,]

  # combine in table format
  qaly_gain_all <- merge(qaly_gain_mean,qaly_gain_CrI_min,by='label',suffixes = c('_mean','_CrI_min'))
  names(qaly_gain_CrI_max)[-1] <- paste0(names(qaly_gain_CrI_max)[-1],'_CrI_max')
  qaly_gain_all <- merge(qaly_gain_all,qaly_gain_CrI_max,by='label',suffixes = c('','_CrI_max'))
  
  # combine in string format "mean [min;max]"
  qaly_gain_summary <- qaly_gain_mean
  for(i in 2:ncol(qaly_gain_mean)){
    qaly_gain_summary[,i] <- paste0(format_table(qaly_gain_mean[,i]),' [',format_table(qaly_gain_CrI_min[,i]),';',format_table(qaly_gain_CrI_max[,i]),']')
  }
  qaly_gain_summary
  qaly_gain_summary$n_reinfect
  qaly_gain_summary$n_infect
  
  # write files
  write.table(qaly_gain_summary,file=paste0(file_tag,'.csv'),sep=',',dec='.',row.names = FALSE)
  
  # summary
  qaly_gain_all <-  merge(qaly_gain_mean,qaly_gain_CrI_min,by='label',suffixes = c('','_CrI_min'))
  qaly_gain_all <-  merge(qaly_gain_all,qaly_gain_CrI_max,by='label') # names(qaly_gain_CrI_max) already updated
  names(qaly_gain_all)
  saveRDS(qaly_gain_all,file=paste0(file_tag,'.rds'))
  
  # return results
  return(list(qaly_gain_all,qaly_gain_summary))  
}

format_table <- function(x){
  return(round(x))
}

get_qaly_loss_daly <- function(output_files,
                          scen_tag = '',
                          time_horizon = NA){
  
  ################################################################ #
  ## SETTINGS AND PREAMBLE ----
  ################################################################ #
  
  # COVID-19 health states and disutility.
  # ref: "Whittington et al. 2022. Value in Health. The Cost-Effectiveness of Remdesivir for Hospitalized Patients With COVID-19"
  disutility_symptomatic       <- 0.19
  disutility_incr_hosp_general <- 0.30 # Disutility of COVID hospitalization with or without supplemental oxygen
  disutility_incr_hosp_oxigen  <- 0.30 # Disutility of COVID hospitalization requiring non-invasive ventilation or high- flow oxygen
  disutility_incr_hosp_mechanical_ventilation <- 0.60    #Disutility of COVID hospitalization requiring mechanical ventilation or ECMO
  disutility_long_covid        <- NA 
  
  # set QALY for SCM health states
  qaly_loss_moderate  <- disutility_symptomatic
  qaly_loss_severe    <- disutility_symptomatic
  qaly_loss_hospital  <- (disutility_symptomatic + disutility_incr_hosp_general)
  qaly_loss_icu       <- (disutility_symptomatic + disutility_incr_hosp_mechanical_ventilation)
  qaly_loss_recovered <- 0
  
  # discounting
  discount_effects <- 0.015
  
  # age-specific QoL
  be_qol <- popnormINDEX(age = 15:100,sex='B',region='BE',year=2018)
  be_qol_age <- c(rep(be_qol$mean[1],15), #0:14
                  be_qol$mean)      #15:100        
  
  # life expectancy
  life_table_be <- get_be_life_table(f_year="2019",f_disc_rate_qaly = discount_effects,qol_age=be_qol_age)
  life_table_be
  
  # use mid-age value: 0:9 ==> 4, 10:19 ==> 14, etc.
  life_table_be[seq(5,100,10),]
  life_table_scm <- life_table_be[seq(5,100,10),]
  life_table_scm
  
  # convert daily into yearly rates
  num_days_year <- 365
  
  ################################################################ #
  ## MODEL OUTPUT ----
  ################################################################ #
  # model output 
  files <- output_files[grepl(scen_tag,output_files)]
  
  if(length(files) == 0 ){
    return(rep(0,9))
  }
  
  prev_I_asy    <- readRDS(files[grepl('prev_I_asymp_sims',files)])
  prev_I_mild   <- readRDS(files[grepl('prev_I_mild_sims',files)])
  prev_I_sev    <- readRDS(files[grepl('prev_I_severe_sims',files)])
  prev_I_hosp   <- readRDS(files[grepl('prev_I_hosp_sims',files)])
  prev_I_icu    <- readRDS(files[grepl('prev_I_icu_sims',files)])
  dim(prev_I_icu)
  
  inc_mortality <- readRDS(files[grepl('mortality_age',files)])
  inc_mortality <- inc_mortality[,,-1] # remove time step column
  dim(inc_mortality)
  
  # model start and dates
  model_date_start <- get_scm_start_date()
  scm_dates        <- model_date_start + (1:ncol(prev_I_asy)-1)
  
  # time horizon
  #time_horizon <- c(get_scm_start_date(),'2022-01-01')
  if(!any(is.na(time_horizon))){
    scm_date_selection <- scm_dates >= as.Date(time_horizon[1]) & scm_dates < as.Date(time_horizon[2])
  } else{
    scm_date_selection <- scm_dates >= model_date_start # always TRUE
  }
  prev_I_asy    <- prev_I_asy[,scm_date_selection,]
  prev_I_mild   <- prev_I_mild[,scm_date_selection,]
  prev_I_sev    <- prev_I_sev[,scm_date_selection,]
  prev_I_hosp   <- prev_I_hosp[,scm_date_selection,]
  prev_I_icu    <- prev_I_icu[,scm_date_selection,]
  dim(prev_I_icu)
  
  inc_mortality <- inc_mortality[,scm_date_selection,] # remove time step column
  dim(inc_mortality)
  
  dim(prev_I_mild)
  dim(prev_I_hosp)
  dim(prev_I_icu) 
  
  # mortality (age-specific)
  dim(inc_mortality)
  sum_mort_age <- apply(inc_mortality,3,rowSums)
  dim(sum_mort_age)
  
  get_lly <- function(life_expectancy_age,sum_mort_age){
    life_expectancy_matrix      <- matrix(rep(life_expectancy_age,nrow(sum_mort_age)),
                                          ncol=ncol(sum_mort_age),
                                          byrow = T)
    lost_life_years             <- sum_mort_age * life_expectancy_matrix
    return(rowSums(lost_life_years))
  } 
  
  lost_life_years         <- get_lly(life_expectancy_age = life_table_scm$life_expectancy,sum_mort_age=sum_mort_age)
  lost_life_years_disc    <- get_lly(life_expectancy_age = life_table_scm$life_expectancy_disc,sum_mort_age=sum_mort_age)
  
  # lost quality-adjusted life years
  lost_qaly       <- get_lly(life_expectancy_age = life_table_scm$life_expectancy_qol,sum_mort_age=sum_mort_age)
  lost_qaly_disc  <- get_lly(life_expectancy_age = life_table_scm$life_expectancy_qol_disc,sum_mort_age=sum_mort_age)
  
  
  # morbidity (not age-specific)
  sum_prev_I_mild <- colSums(apply(prev_I_mild,1,rowSums))
  sum_prev_I_sev  <- colSums(apply(prev_I_sev,1,rowSums))
  sum_prev_I_hosp <- colSums(apply(prev_I_hosp,1,rowSums))
  sum_prev_I_icu  <- colSums(apply(prev_I_icu,1,rowSums))
  sum_inc_mort    <- colSums(apply(inc_mortality,1,rowSums))
  
  # disease related QALY loss
  qaly_loss <- data.table(mild     = (sum_prev_I_mild/num_days_year * qaly_loss_moderate),
                          sev      = (sum_prev_I_sev/num_days_year * qaly_loss_severe),
                          hosp     = (sum_prev_I_hosp/num_days_year * qaly_loss_hospital),
                          icu      = (sum_prev_I_icu/num_days_year * qaly_loss_icu),
                          lly      = lost_life_years,
                          lly_disc = lost_life_years_disc,
                          lqaly    = lost_qaly,
                          lqaly_disc = lost_qaly_disc)
  
  # aggregate
  qaly_loss$total      <- rowSums(qaly_loss[,c(1:4,7)]) # quality-of-life adjusted loss
  qaly_loss$total_disc <- rowSums(qaly_loss[,c(1:4,8)]) # quality-of-life adjusted and discounted loss
  
  # add number of deaths
  qaly_loss$n_mort <- sum_inc_mort
  
  # return results
  return(qaly_loss)
  #print(head(qaly_loss))
}


# based on Sandmann et al, Lancet ID
get_qaly_loss <- function(output_files,
                          scen_tag = '',
                          time_horizon = NA){
  
  ################################################################ #
  ## SETTINGS AND PREAMBLE ----
  ################################################################ #
  
  # COVID-19 health states and QALY LOSS
  # ref: "The potential health and economic value of SARS-CoV-2 vaccination alongside physical distancing in the UK: a transmission model-based future scenario analysis and economic evaluation. Lancet Infect Dis 2021;9"
  # QALY loss per non-fatal case
  qaly_loss_moderate  <- 0.008
  qaly_loss_hospital  <- 0.021
  qaly_loss_icu       <- 0.15
  qaly_loss_recovered <- 0     # long covid?
  
  # discounting
  discount_effects <- 0.015

  # age-specific QoL
  be_qol <- popnormINDEX(age = 15:100,sex='B',region='BE',year=2018)
  be_qol_age <- c(rep(be_qol$mean[1],15), #0:14
                        be_qol$mean)      #15:100        

  # life expectancy
  life_table_be <- get_be_life_table(f_year="2019",f_disc_rate_qaly = discount_effects,qol_age=be_qol_age)
  life_table_be
  
  # use mid-age value: 0:9 ==> 4, 10:19 ==> 14, etc.
  life_table_be[seq(5,100,10),]
  life_table_scm <- life_table_be[seq(5,100,10),]
  life_table_scm
  
  ################################################################ #
  ## MODEL OUTPUT ----
  ################################################################ #
  # model output 
  files <- output_files[grepl(scen_tag,output_files)]
  
  if(length(files) == 0 ){
    warning('scenario tag not found')
    return(rep(0,9))
  }
  if(length(files[grepl('scen_mild_cases_incr',files)])>1){
    warning('NO QALY OUTPUT: SCENARIO TAG "',scen_tag, '" IS NOT UNIQUE!')
    return(rep(0,9))
  }
  
  inc_I_age_all  <- readRDS(files[grepl('scen_cases_incr',files)])[,,-1]
  inc_I_age_mild <- readRDS(files[grepl('scen_mild_cases_incr',files)])[,,-1]
  inc_I_age_hosp <- readRDS(files[grepl('inc_I_hosp_sims',files)])
  inc_I_age_icu  <- readRDS(files[grepl('inc_I_icu_sims',files)])
  
  inc_mortality <- readRDS(files[grepl('mortality_age',files)])
  inc_mortality <- inc_mortality[,,-1] # remove time step column
  dim(inc_mortality)
  
  new_reinfections <- readRDS(files[grepl('new_reinfect',files)])
  dim(new_reinfections)
  
  # model start and dates
  model_date_start <- get_scm_start_date()
  scm_dates        <- model_date_start + (1:ncol(inc_I_age_mild)-1)
  
  # time horizon
  if(!any(is.na(time_horizon))){
    scm_date_selection <- scm_dates >= as.Date(time_horizon[1]) & scm_dates < as.Date(time_horizon[2])
  } else{
    scm_date_selection <- scm_dates >= model_date_start # always TRUE
  }
  inc_I_age_all    <- inc_I_age_all[,scm_date_selection,]
  inc_I_age_mild   <- inc_I_age_mild[,scm_date_selection,]
  inc_I_age_hosp   <- inc_I_age_hosp[,scm_date_selection,]
  inc_I_age_icu    <- inc_I_age_icu[,scm_date_selection,]
  dim(inc_I_age_icu)
  
  inc_mortality <- inc_mortality[,scm_date_selection,] # remove time step column
  dim(inc_mortality)
  
  new_reinfections <- new_reinfections[scm_date_selection,]
  
  # mortality (age-specific)
  sum_mort_age <- apply(inc_mortality,3,rowSums)
  get_lly <- function(life_expectancy_age,sum_mort_age){
    life_expectancy_matrix      <- matrix(rep(life_expectancy_age,nrow(sum_mort_age)),
                                          ncol=ncol(sum_mort_age),
                                          byrow = T)
    lost_life_years             <- sum_mort_age * life_expectancy_matrix
    return(lost_life_years)
  } 
  
  lost_life_years         <- get_lly(life_expectancy_age = life_table_scm$life_expectancy,sum_mort_age=sum_mort_age)
  lost_life_years_disc    <- get_lly(life_expectancy_age = life_table_scm$life_expectancy_disc,sum_mort_age=sum_mort_age)
  
  # lost quality-adjusted life years
  lost_qaly       <- get_lly(life_expectancy_age = life_table_scm$life_expectancy_qol,sum_mort_age=sum_mort_age)
  lost_qaly_disc  <- get_lly(life_expectancy_age = life_table_scm$life_expectancy_qol_disc,sum_mort_age=sum_mort_age)
  
  # morbidity (age-specific)
  sum_inc_I_age_all  <- apply(inc_I_age_all,3,rowSums)
  sum_inc_I_age_mild <- apply(inc_I_age_mild,3,rowSums)
  sum_inc_I_age_hosp <- apply(inc_I_age_hosp,3,rowSums)
  sum_inc_I_age_icu  <- apply(inc_I_age_icu,3,rowSums)
  sum_inc_mort_age   <- sum_mort_age
  
  # morbidity (not age-specific)
  sum_inc_I_all  <- rowSums(sum_inc_I_age_all)
  sum_inc_I_mild <- rowSums(sum_inc_I_age_mild)
  sum_inc_I_hosp <- rowSums(sum_inc_I_age_hosp)
  sum_inc_I_icu  <- rowSums(sum_inc_I_age_icu)
  sum_inc_mort   <- rowSums(sum_inc_mort_age)
  
  # reinfections (not age-specific)
  sum_reinfections <- colSums(new_reinfections)
  
  # get mild infections: mild episodes without hospital admission
  sum_inc_I_mild     <- sum_inc_I_mild - (sum_inc_I_hosp + sum_inc_I_icu)
  sum_inc_I_age_mild <- sum_inc_I_age_mild - (sum_inc_I_age_hosp + sum_inc_I_age_icu)
  
  # disease-related QALY loss
  qaly_loss <- data.table(mild     = (sum_inc_I_mild * qaly_loss_moderate),
                          sev      = 0,
                          hosp     = (sum_inc_I_hosp * qaly_loss_hospital),
                          icu      = (sum_inc_I_icu * qaly_loss_icu),
                          lly      = rowSums(lost_life_years),
                          lly_disc = rowSums(lost_life_years_disc),
                          lqaly    = rowSums(lost_qaly),
                          lqaly_disc = rowSums(lost_qaly_disc),
                          
                          # add number of deaths and cases
                          n_infect  = sum_inc_I_all,
                          n_mild = sum_inc_I_mild,
                          n_hosp = sum_inc_I_hosp,
                          n_icu  = sum_inc_I_icu,
                          n_mort = sum_inc_mort,
                          n_reinfect = sum_reinfections,
                          
                          # add number of deaths and cases
                          n_infect_age  = sum_inc_I_age_all,
                          n_mild_age = sum_inc_I_age_mild,
                          n_hosp_age = sum_inc_I_age_hosp,
                          n_icu_age  = sum_inc_I_age_icu,
                          n_mort_age = sum_inc_mort_age,
                          
                          # age- and severity specific QALY loss
                          qaly_mild_age = sum_inc_I_age_mild * qaly_loss_moderate,
                          qaly_hosp_age = sum_inc_I_age_hosp * qaly_loss_hospital,
                          qaly_icu_age  = sum_inc_I_age_icu * qaly_loss_icu,
                          qaly_mort_age = lost_qaly)
  
  dim(qaly_loss)
  names(qaly_loss) <- gsub('.V','',names(qaly_loss))
  names(qaly_loss)
  
  # aggregate
  qaly_loss$illness    <- rowSums(qaly_loss[,c(1:4)])
  qaly_loss$total      <- qaly_loss$illness + qaly_loss$lqaly # quality-of-life adjusted loss
  qaly_loss$total_disc <- qaly_loss$illness + qaly_loss$lqaly_disc # quality-of-life adjusted and discounted loss
  
  # percentage QALY loss due to mortality
  qaly_loss$perc_mort <- qaly_loss$lqaly / qaly_loss$total * 100
  qaly_loss$perc_mort_disc <- qaly_loss$lqaly_disc / qaly_loss$total * 100
  
  # percentage QALY loss due to mild disease
  qaly_loss$perc_mild <- qaly_loss$mild / qaly_loss$total * 100
  qaly_loss$perc_mild_disc <- qaly_loss$mild / qaly_loss$total * 100
  
  # return results
  return(qaly_loss)
}

# Use data from STATBEL
f_disc_rate_qaly <- 0.03
f_year = 2020
qol_age <- rep(1,100)
get_be_life_table <- function(f_year = 2020 ,f_disc_rate_qaly,qol_age_vector = NA) {
  
  life_table_statbel <- read.xlsx('data/statbel_sterftetafelsAE.xlsx',sheet=paste(f_year),startRow = 3)
  dim(life_table_statbel)
  names(life_table_statbel)
  # WARNING: NOT UNIQUE, some are duplicated for {male, female, both}
  
  # column selection 
  life_table <- life_table_statbel[1:101,c(2,23)]
  names(life_table) <- c('age','life_expectancy')
  tail(life_table)
  
  life_table$life_expectancy_disc     <- get_discounted_life_expectancy(life_expectancy_vector = life_table$life_expectancy, 
                                                                        f_disc_rate_qaly = f_disc_rate_qaly)
  
  life_table$life_expectancy_qol      <- get_discounted_life_expectancy(life_expectancy_vector = life_table$life_expectancy, 
                                                                        f_disc_rate_qaly = 0,
                                                                        qol_age_vector = qol_age_vector)
  
  life_table$life_expectancy_qol_disc <- get_discounted_life_expectancy(life_expectancy_vector = life_table$life_expectancy, 
                                                                        f_disc_rate_qaly = f_disc_rate_qaly,
                                                                        qol_age_vector = qol_age_vector)
  return(life_table)
}

# option to adjust for a given age-specific Quality-of-Life (QoL)
get_discounted_life_expectancy <- function(life_expectancy_vector, 
                                           f_disc_rate_qaly,
                                           qol_age_vector = NA){

  # make sure the QoL vector is set
  if(any(is.na(qol_age_vector))){
    qol_age_vector <- rep(1,length(life_expectancy_vector))
  }
  
  # make sure the QoL vector has the correct dimensions
  if(length(life_expectancy_vector) > length(qol_age_vector)){
    qol_age_vector <- c(qol_age_vector,
                        rep(qol_age_vector[length(qol_age_vector)],
                            length(life_expectancy_vector) - length(qol_age_vector)))
    warning("extended the age-specific Quality-of-Life vector")
  }
  
  # initialize the output vector
  life_expectancy_qol_disc <- vector(length=length(life_expectancy_vector))
  
  # loop over the ages, and apply discounting
  i <- 1
  for(i in 1:length(life_expectancy_vector)){
    
    tmp_years_full      <- 1:floor(life_expectancy_vector[i])
    tmp_year_incomplete <- ceiling(life_expectancy_vector[i])
    tmp_year_remaining  <- life_expectancy_vector[i] - floor(life_expectancy_vector[i])
    
    # QoL-adjusted and discounted life expectancy
    life_expectancy_qol_disc_floor     <- sum(qol_age_vector[(i-1)+tmp_years_full]/((1+f_disc_rate_qaly)^tmp_years_full))
    life_expectancy_qol_disc_remaining <- sum(qol_age_vector[tmp_year_incomplete]/((1+f_disc_rate_qaly)^tmp_year_incomplete)) * tmp_year_remaining
    life_expectancy_qol_disc[i]        <- life_expectancy_qol_disc_floor + life_expectancy_qol_disc_remaining
  }
  
  return(life_expectancy_qol_disc)
}

