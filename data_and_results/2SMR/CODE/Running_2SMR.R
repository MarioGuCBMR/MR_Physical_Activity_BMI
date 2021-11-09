##############
#INTRODUCTION#
##############

#This is a code to produce MR results and QC by loading an already prepared 2SMR dataframe.

###################
#Loading libraries#
###################

library(TwoSampleMR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)
library(rmarkdown)
library(data.table)
library(jsonlite)
library(httr)
library(tidyverse)

###################
#Loading functions#
###################

mr_plots <- function(dat)
{
  #A function to plot sensitivity plots all together.
  
  require(TwoSampleMR)
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  
  temp <- subset(dat, outcome == outcome[1] & exposure == exposure[1])
  exposure_name <- temp$exposure[1]
  outcome_name <- temp$outcome[1]
  
  if(! "labels" %in% names(dat)) dat$labels <- NA
  
  exposure_units <- temp$units.exposure[1]
  outcome_units <- temp$units.outcome[1]
  
  mrs <- mr_singlesnp(temp, all_method=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  mrl <- mr_leaveoneout(temp)
  
  mrs$index <- 1:nrow(mrs)
  mrl$index <- 1:nrow(mrl)
  
  mrs <- dplyr::arrange(merge(mrs, select(temp, SNP, labels), all.x=TRUE), index)
  mrl <- dplyr::arrange(merge(mrl, select(temp, SNP, labels), all.x=TRUE), index)
  
  mrres <- mr(temp, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  
  gridExtra::grid.arrange(
    mr_forest_plot(mrs)[[1]] +
      ggplot2::labs(
        title="a)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")")
      ),
    mr_scatter_plot(mrres, temp)[[1]] +
      ggplot2::labs(
        title="b)",
        x=paste0("SNP effect on ", exposure_name),
        y=paste0("SNP effect on ", outcome_name)
      ) +
      geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    mr_leaveoneout_plot(mrl)[[1]] +
      ggplot2::labs(
        title="c)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")"),
        y="Excluded variant"
      ),
    mr_funnel_plot(mrs)[[1]] +
      ggplot2::labs(title="d)") +
      ggplot2::theme(legend.position="none") +
      ggrepel::geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    ncol=2
  )
}

##############
#Loading data#
##############

path_ <- #write path name finishing with a "/",
file_name <- "" #write the path to your data
file <- paste(path_, file_name, sep = "")

dat_1_pre <- readRDS(file)

###############################################
#Checking mF and Isq and chromose and position#
###############################################

F_ = ((dat_1_pre$beta.exposure)^2)/((dat_1_pre$se.exposure)^2)
mF  = mean(F_)

print(mF)

Isq(dat_1_pre$beta.exposure, dat_1_pre$se.exposure)

####################
#Heterogeneity test#
####################

library(meta)

dat_1_pre$mr <- dat_1_pre$beta.outcome/dat_1_pre$beta.exposure
dat_1_pre$mr_se <- ((dat_1_pre$mr*((dat_1_pre$se.exposure/dat_1_pre$beta.exposure)^2+(dat_1_pre$se.outcome/dat_1_pre$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_pre$mr, dat_1_pre$mr_se)

###################
#First MR analysis#
###################

mr(dat_1_pre)

#######################
#First Rucker Analysis#
#######################

mr_rucker(dat_1_pre)

#########################
#First sensitivity plots#
#########################

dat_1_pre$labels <- NA

image_file <- paste(path_, file_name, ".png", sep = "")

tiff("image_file", units="in", width=10, height=10, res=300)
mr_plots(dat_1_pre)
dev.off()

###################
#REMOVING OUTLIERS#
###################

library(RadialMR)

radial_input <- format_radial(dat_1_pre$beta.exposure, dat_1_pre$beta.outcome, dat_1_pre$se.exposure, dat_1_pre$se.outcome, dat_1_pre$SNP)

radial_output <- egger_radial(radial_input, 0.05, 3) 

outliers <- radial_output$outliers$SNP 

dat_1_post <- dat_1_pre[which(!(dat_1_pre$SNP%in%outliers)),]

##########################################
#Calculating standard QC for the variants#
##########################################

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)
mF  = mean(F_)

print(mF)

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)

##########################
#First heterogeneity test#
##########################

library(meta)

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_post$mr, dat_1_post$mr_se)

#####################
#Generating analysis#
#####################

mr(dat_1_post)

####################
#CALCULATING RUCKER#
####################

mr_rucker(dat_1_post)

#We are gonna double-check:

mr_rucker_jackknife(dat_1_post)

#We checked with a jackknife and it seems that 
#the model for IVW is the best.

#################
#PRINTING IMAGES#
#################

dat_1_post$labels <- NA

image_file <- paste(path_, file_name, "_after_outlier_extraction", ".png", sep = "")

tiff("image_file", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()
