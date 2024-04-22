########################################################################### #
# This file is part of the Stochastic Compartmental Model for SARS-COV-2 
# transmission in Belgium, conceived by members of SIMID group during the 
# COVID19 pandemic.
#
# This file is used to plot the figures for the VoC manuscript
#
# Copyright 2024, SIMID                                        
########################################################################### #

rm(list=ls())

# load functions and data
source('R/main_vaccination.R')
source('R/lib_qaly.R')

################################################################ #
## SETTINGS AND PREAMBLE ----
################################################################ #

# SIMULATION DATA
output_dir <- "output/v28_uptake/"

# select scenario option 
prefix_subdir <- 'orig'
# prefix_subdir <- 'omicron'

# select parameter set
sel_param <- 've00_customQ'
# sel_param <- 've30_customQ'
# sel_param <- 've00_aggrQ'

# set plot options
include_copyright <- FALSE
show_stages     <- FALSE

# to present (aggregated) q-parameters
bool_show_q_param      <- FALSE

# set x-axis limits
x_axis_scm_day      <- sim_date2day(c("2021-09-01","2022-03-01"))
scen_reference_date <- '2021-07-01'  # omit from figure

if(prefix_subdir == 'start'){
  x_axis_scm_day[1] <- 0
  scen_reference_date <- '2020-01-01'
  show_stages <- TRUE
}

# use default y axis?
default_y_axis <- TRUE

# load and select output file
output_files <- dir(output_dir,full.names = T,recursive = T)
file_pattern <- paste(sel_param,prefix_subdir,sep='.*')
output_files <- output_files[grepl(file_pattern,output_files)]
length(output_files)

# check output files
if(length(output_files)==0){
  stop('OUTPUT SELECTION ISSUE IN "',output_dir,'" WITH "',file_pattern,'"')
}

# reference data
be_ref_data      <- read.csv('data/covid19_reference_data_20230526_1030_belgium.csv')

# set figure legend and color library
col_lib <- data.frame(tags = c('xxxx',
                               'orig.*u2ndbooster0_belgium',                       # reference
                                paste0(prefix_subdir,'.*u2ndbooster0_belgium'),    # reported uptake
                                paste0(prefix_subdir,'.*ubooster200_.*_belgium'),  # increased booster uptake
                                paste0(prefix_subdir,'.*infant_belgium'),          # increased infant uptake
                                paste0(prefix_subdir,'.*ubooster60_belgium')),     # reduced booster uptake

                      col  = c('grey','darkred',
                               'gold','darkgreen','darkblue','orange3'),
                      label = c(paste0('Reported (2024-01-01)'),
                                'Estimated dynamics with reported vaccine uptake',
                                'Reported vaccine uptake',
                                'Increased adult booster dose uptake',
                                'Increased 5-11y two doses uptake',
                                'Reduced adult booster dose uptake'),
                      lty = 1,
                      lwd = c(NA,2,2,2,2,2),
                      pch = c(1,NA,NA,NA,NA,NA))
bool_file_sort_decreasing <- T

# y_range_lib
y_range_lib <- rbind(c('/hosp_adm',30,1200),
                     c('hosp_load',0,8000),
                     c('icu_load',0,1200),
                     c('scen_mortality_xtra',0,100))

# set default reference scenario
col_lib$cea_reference <- col_lib$col=='darkred'

# change colors to colorblind friendly set
col_lib$col[2:3] <- 'black'
col_lib$col[4:6] <- brewer.pal(3, "Dark2")#[c(1)]
col_lib$lty[1:6] <- c(1,3,3,1,2,4)

if(prefix_subdir == 'orig'){
  col_lib<- col_lib[-3,]
} else{
  col_lib<- col_lib[-2,]
}

if(prefix_subdir %in% c('omicron')){
  col_lib$cea_reference <- col_lib$col=='black'
}

if(prefix_subdir == 'start'){
  y_range_lib[3:4,3]<- c(1550,220)
}

if(bool_show_q_param){
  col_lib <- col_lib[!is.na(col_lib$pch) | grepl('q-parameters',col_lib$label),]
} else{
  col_lib <- col_lib[!grepl('q-parameters',col_lib$label),]
}

if(prefix_subdir == 'omicron'){
  col_lib$label[-1] <- paste0(col_lib$label[-1],' [no Omicron VOC]')
}

color_indices <- get_color_index_multiple(output_files,col_lib = col_lib,show_warnings = FALSE)
table(color_indices)
output_files  <- output_files[!is.na(color_indices)]

sel_col_indices <- unique(color_indices,na.last=T)
sel_col_indices <- sel_col_indices[!is.na(sel_col_indices)]
col_lib         <- col_lib[sort(c(1,sel_col_indices)),]

## MODEL PARAMETERS ----

# load model parameters + wave and scenario info
param_file_names         <- dir(output_dir,pattern = 'MCMCmulti.*',full.names = T,recursive = T)
scm_num_waves            <- vector(length = length(param_file_names))
scm_callibration_day     <- vector(length = length(param_file_names))
scm_region               <- vector(length = length(param_file_names))
scm_scenario_cp_day      <- c()

for(i_file  in 1:length(param_file_names)){
  parms                         <- read.csv(param_file_names[i_file],header=T)
  scm_num_waves[i_file]         <- identify_total_nb_stages(names(parms))
  scm_callibration_day[i_file]  <- unique(parms$ndays_calibration)
  scm_region[i_file]            <- unique(parms$region_id)
  scm_scenario_cp_day           <- c(scm_scenario_cp_day,unique(parms[grepl('cnt_adjust_day',names(parms))]))
}

# number stages in which social contact behaviour changed
scm_num_cp <- unique(scm_num_waves)
if(!show_stages){
  scm_num_cp <- 0
}

# set calibration and scenario dates
scm_callibration_date    <- unique(sim_day2date(scm_callibration_day))
scm_scenario_cp          <- sort(unique(sim_day2date(unlist(scm_scenario_cp_day))))


if(bool_show_q_param){
  opt_aggr_waves  <- list(1:8,9:11,12:17,18:31,32:42)
  scm_scenario_cp_day <- get_CoMix_change_day(unlist(lapply(opt_aggr_waves,min)))
  scm_scenario_cp     <- sort(unique(sim_day2date(unlist(scm_scenario_cp_day))))
  
}

# sort output files
output_files <- sort(output_files,decreasing = bool_file_sort_decreasing)

# REFERENCE DATA ----
be_ref_data      <- read.csv('data/covid19_reference_data_20230526_1030_belgium.csv')
be_ref_data$date <- as.Date(be_ref_data$date)

# contact data
CoMix_mat <- readRDS(file='data/social_contact_data.rds')

## PDF STREAM ----
output_dir_base <- unlist(strsplit(output_dir,'/'))
output_dir_base <- output_dir_base[length(output_dir_base)]

# set output directory
output_dir_summary <- gsub(basename(output_dir),paste(basename(output_dir),'summary',sep='_'),output_dir)
if(!dir.exists(output_dir_summary)) {dir.create(output_dir_summary)}
print(output_dir_summary)

# set file name based on settings above
pdf_name <- paste0('incidence_',prefix_subdir,
                   '_',sel_param,
                   ifelse(min(x_axis_scm_day)==0,'_full',''),
                   ifelse(bool_show_q_param,'_q',''),'.pdf')

print(pdf_name)
pdf(file.path(output_dir_summary,pdf_name),7,4)
par(mar=c(4.5,5,1,1))
x_axis_vaccine = c(300,max(x_axis_scm_day))

# ## POPULATION LEVEL (lines and polygon) ----
multi_plot_incidence_time(output_files,
                          col_lib        = col_lib,
                          bool_polygon_multi = T,
                          x_axis_scm_day = x_axis_scm_day,
                          be_ref_data    = be_ref_data,
                          db_C_sim     = CoMix_mat$db_C_sim,
                          default_y_axis = default_y_axis,
                          y_range_lib    = y_range_lib,
                          include_date_in_legend = FALSE,
                          include_copyright = include_copyright,
                          scm_callibration_date  = NA,
                          scm_scenario_cp        = scm_scenario_cp,
                          scm_num_cp             = scm_num_cp)


# # BY AGE: bars and polygon (time consuming) ----
# multi_plot_incidence_age_time(output_files = output_files,
#                               col_lib = col_lib,
#                               x_axis_scm_day = x_axis_scm_day,
#                               x_axis_vaccine = x_axis_vaccine)
dev.off()

pdf(file.path(output_dir_summary,gsub('.pdf','_single.pdf',pdf_name)),7,4)#,9,4)
par(mar=c(4.5,5,1,1))

multi_plot_incidence_time(output_files[grepl('/hosp_adm',output_files)],
                          col_lib        = col_lib,
                          bool_polygon_multi = T,
                          x_axis_scm_day = c(sim_date2day('2021-10-15'),x_axis_scm_day[2]),
                          be_ref_data    = be_ref_data,
                          db_C_sim       = CoMix_mat$db_C_sim,
                          default_y_axis = FALSE,
                          y_range_lib    = y_range_lib,
                          include_date_in_legend = FALSE,
                          include_copyright = include_copyright,
                          bool_x_lable_full = TRUE,
                          scm_callibration_date  = NA,
                          scm_scenario_cp        = scm_scenario_cp,
                          scm_num_cp             = scm_num_cp)

add_vertical_line('2021-12-25',bool_text =T,date_tag= 'Dominance Omicron',pos_factor=1.0)
dev.off()  


# # CUMULATIVE ----
pdf(file.path(output_dir_summary,gsub('incidence_','cumulative_',pdf_name)),7,4)#,9,4)
par(mar=c(4.5,5,1,1))
multi_plot_cumulative_time(output_files,
                          col_lib               = col_lib, 
                          x_axis_scm_day        = x_axis_scm_day,
                          scen_reference_date   = scen_reference_date,
                          db_C_sim              = CoMix_mat$db_C_sim,
                          scm_callibration_date = scm_callibration_date,
                          scm_scenario_cp       = scm_scenario_cp,
                          scm_num_cp            = scm_num_cp)
dev.off()

## QALY ----

if(any(col_lib$cea_reference)){
  
  pdf_name_qaly <- file.path(output_dir_summary,gsub('incidence_','qaly_gain_',pdf_name))
  pdf(pdf_name_qaly,10,4)
  
  multi_qaly_loss(output_files =  output_files ,
                  col_lib = col_lib,
                  i_row_reference = which(col_lib$cea_reference),
                  time_horizon = range(sim_day2date(c(0,x_axis_scm_day))),
                  file_tag  = gsub('.pdf','',pdf_name_qaly),
                  bool_plot = TRUE) -> qaly_out
  dev.off()
}

# COMBINED QALY FIGURE ####
# note: this is not executed by default 
if(0==1){ 
  output_dir_summary <- "output/v28_uptake_summary"

  qaly_gain_orig    <- readRDS(dir(output_dir_summary,pattern = 'qaly_gain_orig_ve00_customQ.*rds',full.names = T))
  qaly_gain_omicron <- readRDS(dir(output_dir_summary,pattern = 'qaly_gain_omicron_ve00_customQ.*rds',full.names = T))
  qaly_gain_ve30    <- readRDS(dir(output_dir_summary,pattern = 'qaly_gain_orig_ve30_customQ.*rds',full.names = T))
  
  col_names <- c('label','total','total_CrI_min','total_CrI_max')
  qaly_gain_orig    <- qaly_gain_orig[,col_names]
  qaly_gain_omicron <- qaly_gain_omicron[,col_names]
  qaly_gain_ve30    <- qaly_gain_ve30[,col_names]
  
  qaly_gain_omicron$label <- gsub('\\[.*','',qaly_gain_omicron$label)
  qaly_gain_omicron$label <- gsub('Reported vaccine uptake','Estimated dynamics with reported vaccine uptake',qaly_gain_omicron$label)

  qaly_gain_orig    <- qaly_gain_orig[order(qaly_gain_orig$label),]
  qaly_gain_ve30    <- qaly_gain_ve30[order(qaly_gain_ve30$label),]
  qaly_gain_omicron <- qaly_gain_omicron[order(qaly_gain_omicron$label),]
  
   add_qaly_arrows <- function(qaly_gain_matrix,y_sel,y_dev,qaly_col){
     
    y_ticks <- 1:length(y_sel)+y_dev
    arrows(x0=qaly_gain_matrix$total_CrI_min[y_sel],
           y0=y_ticks,
           x1=qaly_gain_matrix$total_CrI_max[y_sel],
           y1=y_ticks,
           col=qaly_col,
           lwd=2,
           angle = 90,
           length = 0.05)
    arrows(x1=qaly_gain_matrix$total_CrI_min[y_sel],
           y1=y_ticks,
           x0=qaly_gain_matrix$total_CrI_max[y_sel],
           y0=y_ticks,
           col=qaly_col,
           lwd=2,
           angle = 90,
           length = 0.05)
    points(qaly_gain_matrix$total[y_sel],
           y = y_ticks,
           col = qaly_col,
           lwd=2,
           pch = 4)
  }
  
   y_sel <- c(4,3,5)
   pdf(file.path(output_dir_summary,'qaly_gain_wiskers.pdf'),7,4)#,9,4)
   par(mar=c(5,15,1,1))
   plot(x=range(qaly_gain_orig[-1,-1],qaly_gain_ve30[-1,-1]),
        y = range(y_sel+1,y_sel-1),
        ylim = c(length(y_sel)+0.5,0.5),
        col=0,
        yaxt='n',
        ylab='',
        xlab='QALY gain (mean and 95% CrI)')
   
   y_labels = qaly_gain_orig$label[y_sel]
   axis(2,at=1:length(y_labels),labels = y_labels,las=2,lwd=0)
   
    grid(nx=NULL,ny=NA)
    grid(nx=NA,ny=length(y_sel))
    abline(v=0,lty=3)
    add_qaly_arrows(qaly_gain_orig,y_sel,y_dev=-0.2,qaly_col = 'yellow3')
    add_qaly_arrows(qaly_gain_ve30,y_sel,y_dev=0,qaly_col = 'lightblue')
    add_qaly_arrows(qaly_gain_omicron,y_sel,y_dev=+0.2,qaly_col = 'darkgreen')
  
    legend('bottomleft',
           c('VE transmission = 0 and Omicron VOC',
             'VE transmission = 30%',
             'No Omicron VOC'),
           col= c('yellow3','lightblue','darkgreen'),
           lwd=2,
           cex=0.75,
           xpd=TRUE,
           inset=c(-0.75,-0.2)) 
    dev.off()
} # end if-clause for combined QALY figure


# UPTAKE FIGURE ####
# note: this is not executed by default 
if(0==1){  
  
  # load data
  uptake_orig  <- read.csv('data/uptake/uptake_vac_voc/vaccine_uptake_booster_vSCENjan09_u2ndbooster0_belgium.csv')
  uptake_adult <- read.csv('data/uptake/uptake_vac_voc/vaccine_uptake_booster_vSCENjan09_u2ndbooster0_ubooster200_jun01_belgium.csv')
  uptake_child <- read.csv('data/uptake/uptake_vac_voc/vaccine_uptake_booster_vSCENjan09_u2ndbooster0_INFANT_belgium.csv')
  uptake_red   <- read.csv('data/uptake/uptake_vac_voc/vaccine_uptake_booster_vSCENjan09_u2ndbooster0_ubooster60_belgium.csv')
  
  uptake_orig  <- aggregate(. ~ date, data = uptake_orig, sum)
  uptake_adult <- aggregate(. ~ date, data = uptake_adult, sum)
  uptake_child <- aggregate(. ~ date, data = uptake_child, sum)
  uptake_red   <- aggregate(. ~ date, data = uptake_red, sum)
  
  col_1dose         <- grepl('_A_',names(uptake_child)) 
  col_2dose         <- grepl('_B_',names(uptake_child)) 
  col_child_0       <- grepl('_0_',names(uptake_child))
  col_child_10      <- grepl('_10_19',names(uptake_child))
  col_child         <- col_child_0 | col_child_10
  col_booster       <- grepl('_E_',names(uptake_child))
  names(uptake_child)[col_2dose]
  
  # plot
  x_ticks <- pretty(c(as.Date("2021-07-01"),sim_day2date(x_axis_scm_day[2])),10);x_ticks
  x_ticks <- pretty(as.Date(c("2021-08-01",'2022-03-01')),10);x_ticks
  y_scale <- 1e6
  sel_pch <- seq(1,nrow(uptake_child),20)
  
  point_selection <- sim_date2day(add_uptake_boxes(bool_return= TRUE)) - sim_date2day(min(uptake_orig$date))+1
  plot_lwd <- 2
  
  pdf(file.path(output_dir_summary,'scenario_uptake.pdf'),9,5)
  par(mar=c(4.6,5,1,3),mfrow=c(2,1))
  
  # ADULT BOOSTER
  plot(as.Date(uptake_orig$date),
       cumsum(rowSums(uptake_orig[col_2dose & !col_child]))/y_scale,
       type='l',
       xlim=range(x_ticks),
       col='white',
       lwd=2,cex.axis=0.9,
       ylim=c(-5e5,10.5e6)/y_scale,
       xaxt='n',las=2,cex.lab=0.8,
       yaxt='n',
       xlab='',
       ylab='Uptake +20y\n (million doses)')
  grid(nx=NA,ny=NULL)
  axis(1,at=x_ticks,labels=format(x_ticks,'%d/%m/%y'),cex.axis=0.8)
  
  be_pop      <- get_regional_pop()
  y_ticks     <- pretty(c(0,10e6)/y_scale)
  axis(2,at=y_ticks,labels=y_ticks,las=2,cex.axis=0.8)
  
  y_ticks_rel <- y_ticks * y_scale / sum(be_pop[-(1:2)])
  y_ticks     <- y_ticks[y_ticks_rel<1]
  y_ticks_rel <- y_ticks_rel[y_ticks_rel<1]
  axis(4,at=y_ticks,labels=round(y_ticks_rel,2),las=2,cex.axis=0.8)
  axis(4,at=10,labels='%',las=2,cex.axis=0.8,lty = 0)
  
  add_line_with_marks(x_values = uptake_orig$date,
                      y_values = cumsum(rowSums(uptake_orig[col_2dose & !col_child])),
                      point_selection = point_selection,
                      plot_color = col_lib$col[2],
                      plot_pch = col_lib$lty[2],
                      y_scale = y_scale,
                      plot_lwd = plot_lwd)
  add_line_with_marks(x_values = uptake_adult$date,
                      y_values = cumsum(rowSums(uptake_adult[col_booster & !col_child])),
                      point_selection = point_selection,
                      plot_color = col_lib$col[3],
                      plot_pch = col_lib$lty[3],
                      y_scale = y_scale,
                      plot_lwd = plot_lwd)
  add_line_with_marks(x_values = uptake_red$date,
                      y_values = cumsum(rowSums(uptake_red[col_booster & !col_child])),
                      point_selection = point_selection,
                      plot_color = col_lib$col[5],
                      plot_pch = col_lib$lty[5],
                      y_scale = y_scale,
                      plot_lwd = plot_lwd)
  add_line_with_marks(x_values = uptake_orig$date,
                      y_values = cumsum(rowSums(uptake_orig[col_booster & !col_child])),
                      point_selection = point_selection,
                      plot_color = 5,
                      plot_pch = 5,
                      y_scale = y_scale,
                      plot_lwd = plot_lwd)
  
  add_uptake_boxes(y_values=c(0,0,1,1)+9.5,bool_abline = TRUE,opt_sel = 2)
  
  uptake_legend <- data.frame(name=c('2e dose: reported',
                                     '1e booster dose: reported',
                                     '1e booster dose: increased uptake',
                                     '1e booster dose: 60% of reported'),
                              col=c(col_lib$col[2],5,col_lib$col[3],col_lib$col[5]),
                              lty=1,#c(1,1,2,2,2),
                              pch=c(col_lib$lty[2],5,col_lib$lty[3],col_lib$lty[5])
  )
  
  legend('topleft',
         uptake_legend$name,
         col=uptake_legend$col,
         lty=NA,#uptake_legend$lty,
         pch=uptake_legend$pch,
         bg='white',
         #xpd = TRUE,
         adj=0,
         #inset=-0.65,
         ncol=2,
         cex=0.5,
         lwd=plot_lwd)
  
  plot(as.Date(uptake_child$date),
       cumsum(rowSums(uptake_child[col_2dose & col_child]))/y_scale,
       type='l',
       xlim=range(x_ticks),
       col='white',
       lwd=2,cex.axis=0.9,cex.lab=0.8,
       ylim=c(-0.1e6,2.2e6)/y_scale,
       xaxt='n',las=2,
       yaxt='n',
       xlab='Time',
       ylab='Uptake 0-19y\n(million doses)')
  grid(nx=NA,ny=NULL)
  axis(1,at=x_ticks,labels=format(x_ticks,'%d/%m/%y'),cex.axis=0.8)
  
  y_ticks     <- pretty(c(0,2e6)/y_scale)
  axis(2,at=y_ticks,labels=y_ticks,las=2,cex.axis=0.8)
  y_ticks_rel <- y_ticks * y_scale / sum(be_pop[(1:2)])
  y_ticks     <- y_ticks[y_ticks_rel<0.7]
  y_ticks_rel <- y_ticks_rel[1:length(y_ticks)]
  axis(4,at=y_ticks,labels=round(y_ticks_rel,1),las=2,cex.axis=0.8)
  axis(4,at=2,labels='%',las=2,cex.axis=0.8,lty = 0)
  
  add_line_with_marks(x_values = uptake_child$date,
                      y_values = cumsum(rowSums(uptake_child[col_2dose & col_child])),
                      point_selection = point_selection,
                      plot_color = col_lib$col[4],
                      plot_pch = col_lib$lty[4],
                      y_scale = y_scale,
                      plot_lwd = plot_lwd)
  add_line_with_marks(x_values = uptake_orig$date,
                      y_values = cumsum(rowSums(uptake_orig[col_2dose & col_child])),
                      point_selection = point_selection,
                      plot_color = col_lib$col[2],
                      plot_pch = col_lib$lty[2],
                      y_scale = y_scale,
                      plot_lwd = plot_lwd)
  
  add_uptake_boxes(y_values=c(0,0,0.2,0.2)+2,bool_abline = TRUE,opt_sel = 1)
  
  uptake_legend_child <- data.frame(name=c('2e dose: increased',
                                           '2e dose: reported'),
                                    col=col_lib$col[c(4,2)],
                                    pch=col_lib$lty[c(4,2)]
  )
  legend('topright',
         uptake_legend_child$name,
         col=uptake_legend_child$col,
         lty=NA,#uptake_legend$lty,
         pch=uptake_legend_child$pch,
         bg='white',
         adj=0,
         cex=0.5,
         lwd=plot_lwd)
  
  dev.off()
  
}
