########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file contains functions to pre-process the social contact data and to
# define the social contact scenarios.
#
# Copyright 2024, SIMID                                        
########################################################################### #

## HELP FUNCTIONS FOR CONTACT MATRICES

# stage <- 1;bool_asymptomatic = TRUE
select_contact_matrix <- function(contact_data,stage,bool_asymptomatic = TRUE){
  
  # check
  if(stage <1 || stage > nrow(contact_data$db_C_sim)){
    stop("Invalid 'stage' provided to get_contact_matrix")
  }
  
  # get reference matrix
  C_reference <- contact_data$db_C_sim$c_tag[stage]
  
  # specify (a)asymptomatic contact behaviour
  C_reference <- paste(C_reference,ifelse(bool_asymptomatic,'asy','sy'),sep='_')

  return(contact_data[[C_reference]])
}  
  
get_beta_index <- function(beta_list_all, sim_day){
  
  i_beta <- which((sim_day >= beta_list_all$beta_estim_cp_start) &
                     (sim_day <= beta_list_all$beta_estim_cp_end))
  
  if(is.null(i_beta) || length(i_beta)==0){
    return(NULL)
  }
  
  return(i_beta)
}


################################################################ #
## HELP FUNCTIONS W.R.T SOCIAL CONTACTS ----
################################################################ #

# get susceptibility matrices
# note: a delay of 5 means the original behaviour on t=0, starting with a transition 
# from t+1, and full effect by t+5. For example: seq(0,1,length=6)
calculate_all_beta_matrices <- function(contact_matrices, 
                                        q_stage,
                                        behavioral_change_delay,
                                        db_cnt_adjust){  # scenarios
  
  # default
  sel_C_sim            <- contact_matrices$db_C_sim
  sel_C_sim$sim_day    <- sim_date2day(sel_C_sim$date)
  sel_C_sim$q_stage_id <-  as.numeric(sel_C_sim$q_stage_id)
                                       
  # selection
  sel_cnt_adjust <- NULL
  if(!is.null(db_cnt_adjust) && sum(db_cnt_adjust$value)>0){
    sel_cnt_adjust           <- db_cnt_adjust[db_cnt_adjust$value!=0,]
    sel_cnt_adjust$sim_day   <- sim_date2day(sel_cnt_adjust$date)

    if(any(sel_cnt_adjust$is_parmanent)){
      sel_C_sim        <- sel_C_sim[sel_C_sim$sim_day < min(sel_cnt_adjust$sim_day[sel_cnt_adjust$is_parmanent]),]
    }
    
    if(nrow(sel_cnt_adjust)>0){
      # define the corresponding contact data
      sel_cnt_adjust$i_reference      <- NA
      sel_cnt_adjust$reference_day    <- 
        for(i_adjust in 1:nrow(sel_cnt_adjust)){
          i_day <- sim_date2day(sel_cnt_adjust$reference_date[i_adjust])
          sel_cnt_adjust$i_reference[i_adjust] <- max(which(i_day>=sel_C_sim$sim_day))
        }
    }
  } 
  
  # aggregate change points for estimated matrices
  beta_estim_cp_start   <- c(sel_C_sim$sim_day,sel_cnt_adjust$sim_day)
  beta_estim_cp_end     <- beta_estim_cp_start + behavioral_change_delay

  # initiate arrays
  n_stages_total  <- length(beta_estim_cp_start)
  n_stages_data   <- nrow(sel_C_sim)
  beta_estim_asy  <- array(dim = c(10,10,n_stages_total))
  beta_estim_sy   <- array(dim = c(10,10,n_stages_total))
  
  for(i_stage in 1:n_stages_data){
    beta_estim_asy[,,i_stage] <- select_contact_matrix(contact_matrices,i_stage,bool_asymptomatic = TRUE) * q_stage[sel_C_sim$q_stage_id[i_stage],]
    beta_estim_sy[,,i_stage]  <- select_contact_matrix(contact_matrices,i_stage,bool_asymptomatic = FALSE) * q_stage[sel_C_sim$q_stage_id[i_stage],]
  }
  
  if(n_stages_total>n_stages_data){
    for(i_adjust in 1:nrow(sel_cnt_adjust)){
      i_beta_reference          <- sel_cnt_adjust$i_reference[i_adjust]
      i_stage                   <- which(sel_cnt_adjust$sim_day[i_adjust] == beta_estim_cp_start)
      beta_estim_asy[,,i_stage] <- beta_estim_asy[,,i_beta_reference] * sel_cnt_adjust$value[i_adjust]
      beta_estim_sy[,,i_stage]  <- beta_estim_sy[,,i_beta_reference] * sel_cnt_adjust$value[i_adjust]
    }
  }
  
  return(list(beta_estim_asy      = beta_estim_asy,
              beta_estim_sy       = beta_estim_sy,
              beta_estim_cp_start = beta_estim_cp_start,
              beta_estim_cp_end   = beta_estim_cp_end
              ))
}

select_beta <- function(calendar_time,
                        beta_list_all,
                        db_cnt_adjust){
  
  i_beta <- get_beta_index(beta_list_all,calendar_time)
  
  if(is.null(i_beta)){
        return(list(bool_update = FALSE,asy=NA,sy=NA))
  } else
   
     # 1st stage
     if(any(i_beta == 1)){
       beta_asy <- beta_list_all$beta_estim_asy[,,i_beta] 
       beta_sy  <- beta_list_all$beta_estim_sy[,,i_beta] 
     } else if(any(i_beta > 1)){
      
      reference_stage   <- max(i_beta)
      days_diff   <- beta_list_all$beta_estim_cp_end[reference_stage] - beta_list_all$beta_estim_cp_start[reference_stage]
      wStage_prev <- 1 - (calendar_time - beta_list_all$beta_estim_cp_start[reference_stage])/days_diff
      wStage_new  <- 1 - wStage_prev;
      
      beta_asy       <- wStage_prev*beta_list_all$beta_estim_asy[,,reference_stage-1] + wStage_new*beta_list_all$beta_estim_asy[,,reference_stage]
      beta_sy        <- wStage_prev*beta_list_all$beta_estim_sy[,,reference_stage-1] + wStage_new*beta_list_all$beta_estim_sy[,,reference_stage]
    } 
    
  # include safety checks
  if(length(beta_asy)==0){
    stop(paste('STOP: beta_asy is NULL on day', calendar_time))
  }
  if(length(beta_sy)==0){
    stop(paste('STOP: beta_sy is NULL on day', calendar_time))
  }
  if(all(beta_sy == beta_asy)){
    stop(paste('STOP: beta_sy == beta_asy on day', calendar_time))
  }
  
  return(list(bool_update = TRUE,asy=beta_asy,sy=beta_sy))
}



