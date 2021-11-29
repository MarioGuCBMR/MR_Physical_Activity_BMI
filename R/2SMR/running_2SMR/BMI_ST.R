##############
#INTRODUCTION#
##############

#This is a code to run 2SMR with BMI as exposure and sedentary time as outcome.

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

ld_proxies <- function(snp){
  #Enter vector of SNPs it will output all its LD friends from 1KG European LD friends with r2 > 0.8
  #in 500kb.
  
  fake_df <- t(as.data.frame(c("d_prime", "variation2", "population_name", "r2", "variation1")))
  
  colnames(fake_df) <- fake_df[1,]
  rownames(fake_df) <- c(1)
  
  #Setting the server:
  
  server <- "http://grch37.rest.ensembl.org"
  
  for(i in snp){
    
    ext_1 <- paste("/ld/human/", i, sep = "")
    ext_2 <- paste(ext_1, "/1000GENOMES:phase_3:EUR", sep = "")
    
    r <- GET(paste(server, ext_2, sep = ""), content_type("application/json"))
    new <- fromJSON(toJSON(content(r)))
    
    fake_df <- rbind(fake_df, new)
    
  }
  
  #Now filtering for those that are in high LD:
  
  final_df <- fake_df[which(as.numeric(fake_df$r2) > 0.8),] #The NAs by coercion are the rows from the fake_df, ignore them!
  
  return(final_df)
  
}

Isq = function(y,s){
  #A function to calculate I-squared parameter.
  
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#######################
#Loading original data#
#######################

#We are loading the data that we have curated beforehand:

sedentary <- fread("~/CURATED_DATA/Sedentary_doherty.txt")

BMI <- fread("~/CURATED_DATA/BMI_Locke_Chr_Pos.txt")

##############################################
#Getting genome-wide significant SNPs for BMI#
##############################################

BMI_gw <- BMI[which(as.numeric(BMI$p) < 0.00000005),]

summary(BMI_gw$Freq1.Hapmap) #good. #the NAs will be removed. Those above 0.99 too
summary(BMI_gw$p) #good.

#Getting them independent only with OLDrsid

BMI_gw_ind_oldRSID <- BMI_gw

BMI_gw_ind_oldRSID$rsid <- BMI_gw_ind_oldRSID$SNP
BMI_gw_ind_oldRSID$pval <- BMI_gw_ind_oldRSID$p

#Let's check the chromosome first...

check <- BMI_gw_ind_oldRSID[which(stringr::str_detect(BMI_gw_ind_oldRSID$chr_pos, "chr6") == TRUE),]

check <- check[order(check$chr_pos),]

head(check, 30) #those close to the MHC region are above the 34.000.000 borderline.
tail(check, 40) #none.

#Everything is OK.

###########################################
#Comparing whether we have merged variants#
###########################################

which(BMI_gw_ind_oldRSID$SNP != BMI_gw_ind_oldRSID$newRSID) #we have one different!

#We will have to check how does that work with the ld_pruning...

####################################
#LD clumping with the original RSID#
####################################

BMI_gw_ind_oldRSID_ <- ieugwasr::ld_clump_local(BMI_gw_ind_oldRSID, bfile = "~/1K/1kg.v3/1kg.v3/EUR", plink_bin ="~/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We got 68 SNPs, let's check with the new:

BMI_gw_ind_newRSID <- BMI_gw_ind_oldRSID

BMI_gw_ind_newRSID$rsid <- BMI_gw_ind_newRSID$newRSID
BMI_gw_ind_newRSID$pval <- BMI_gw_ind_newRSID$pval

BMI_gw_ind_newRSID_ <- ieugwasr::ld_clump_local(BMI_gw_ind_newRSID, bfile = "~/1K/1kg.v3/1kg.v3/EUR", plink_bin ="~/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

which(!(BMI_gw_ind_newRSID_$chr_pos%in%BMI_gw_ind_oldRSID_$chr_pos)) 

#It seems that the GIANT summary statistics had a merged SNP that is not found on the 1000G reference panel.
#However, if we use the newest SNP, we get one more hit!
#That means that the matching will be done with chromosome and position, but that we will ultimately use the new RSIDs for 
#the analysis.

###########################
#Matching with the outcome#
###########################

sedentary_oldrsid_2_BMI <- sedentary[which(sedentary$chr_pos%in%BMI_gw_ind_oldRSID_$chr_pos),] #68/68 as the classic 2SMR versions.
sedentary_newrsid_2_BMI <- sedentary[which(sedentary$chr_pos%in%BMI_gw_ind_newRSID_$chr_pos),] #68/68 but one SNP has changed.

#The lead SNP changes taking into account which version we are using. 
#We will have to compare, though of course, I doubt the results change a lot since they are bound to be in the same loci.

################################
#Preparing data for TwoSampleMR#
################################

#We are going to do this twice to check whether we wind up with any differences or not:

sedentary_old_set <- sedentary_oldrsid_2_BMI
BMI_old_set <- BMI_gw_ind_oldRSID_

sedentary_old_set$id <- "sedentary"
BMI_old_set$id <- "BMI"

##############################
#FIRST WITH THE ORIGINAL DATA#
##############################

#We should have the same SNPs in each.
#Perfect, now we get the columns for TwoSampleMR

sedentary_4_merge_old <- sedentary_old_set %>%
  select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF, id)

colnames(sedentary_4_merge_old) <- c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "id.outcome")

BMI_4_merge_old <- BMI_old_set %>%
  select(SNP, A1, A2, Freq1.Hapmap, b, se, p, id, N)

colnames(BMI_4_merge_old) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "id.exposure", "samplesize.exposure")

#######################
#Finally we merge them#
#######################

BMI_4_merge_old$exposure <- "BMI"
sedentary_4_merge_old$outcome <- "Sedentary"

dat_1_pre_old <- harmonise_data(BMI_4_merge_old, sedentary_4_merge_old, action = 3)

dat_1_pre_old <- dat_1_pre_old[which(dat_1_pre_old$palindromic == FALSE),]
dat_1_pre_old <- dat_1_pre_old[which(dat_1_pre_old$remove == FALSE),]

#########################################
#WE REPEAT THIS BUT WITH NEWRSID VERSION#
#########################################

#We are going to do this twice to check whether we wind up with any differences or not:

sedentary_new_set <- sedentary_newrsid_2_BMI
BMI_new_set <- BMI_gw_ind_newRSID_

sedentary_new_set$id <- "sedentary"
BMI_new_set$id <- "BMI"

#Perfect, now we get the columns for TwoSampleMR

sedentary_4_merge_new <- sedentary_new_set %>%
  select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF, id)

colnames(sedentary_4_merge_new) <- c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "id.outcome")

BMI_4_merge_new <- BMI_new_set %>%
  select(newRSID, A1, A2, Freq1.Hapmap, b, se, p, id, N)

colnames(BMI_4_merge_new) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "id.exposure", "samplesize.exposure")

#######################
#Finally we merge them#
#######################

BMI_4_merge_new$exposure <- "BMI"
sedentary_4_merge_new$outcome <- "Sedentary"

dat_1_pre_new <- harmonise_data(BMI_4_merge_new, sedentary_4_merge_new, action = 3)

dat_1_pre_new <- dat_1_pre_new[which(dat_1_pre_new$palindromic == FALSE),]
dat_1_pre_new <- dat_1_pre_new[which(dat_1_pre_new$remove == FALSE),]

#We get one more SNP than before.
#Hence, we stick with the new one:

dat_1_pre <- dat_1_pre_new

saveRDS(dat_1_pre, "BMI_ST_MR_DF")

######################
#Let's check the data#
######################

summary(dat_1_pre$eaf.exposure) #perfect
summary(dat_1_pre$eaf.outcome) #perfect.

###############################
#Now let's add the sample size#
###############################

dat_1_pre$samplesize.outcome <- 91105
dat_1_pre$outcome <- "sedentary"
dat_1_pre$exposure <- "BMI"

#Now we can perform steiger:

dat_1_pre <- steiger_filtering(dat_1_pre)

dat_1_pre <- dat_1_pre[which(dat_1_pre$steiger_dir == TRUE),] 

##############################
#THUS we can run our analysis#
##############################

saveRDS(dat_1_pre, "BMI_ST_MR_DF")

###############################################
#Checking mF and Isq and chromose and position#
###############################################

#No need to worry about chromosome 6 here.
#We already checked above and all is good.

F_ = ((dat_1_pre$beta.exposure)^2)/((dat_1_pre$se.exposure)^2)
mF  = mean(F_)

print(mF)

#58.47284

#Cannot find the code for Isq, but we will find it.

Isq(dat_1_pre$beta.exposure, dat_1_pre$se.exposure)

#0.9831416

##########################
#First heterogeneity test#
##########################

library(meta)

dat_1_pre$mr <- dat_1_pre$beta.outcome/dat_1_pre$beta.exposure
dat_1_pre$mr_se <- ((dat_1_pre$mr*((dat_1_pre$se.exposure/dat_1_pre$beta.exposure)^2+(dat_1_pre$se.outcome/dat_1_pre$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_pre$mr, dat_1_pre$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0204 [0.0043; 0.0552]; tau = 0.1430 [0.0656; 0.2349];
#I^2 = 34.4% [10.8%; 51.7%]; H = 1.23 [1.06; 1.44]

#Test of heterogeneity:
#  Q d.f. p-value
#97.50   64  0.0044

###################
#First MR analysis#
###################

mr(dat_1_pre)

#id.exposure id.outcome   outcome exposure                    method nsnp           b         se       pval
#1         BMI  sedentary sedentary      BMI                  MR Egger   65 -0.09695645 0.08848203 0.27734825
#2         BMI  sedentary sedentary      BMI           Weighted median   65  0.03123417 0.04088082 0.44484936
#3         BMI  sedentary sedentary      BMI Inverse variance weighted   65  0.07313790 0.03096465 0.01817759
#4         BMI  sedentary sedentary      BMI               Simple mode   65 -0.05150509 0.10730552 0.63287415
#5         BMI  sedentary sedentary      BMI             Weighted mode   65 -0.05150509 0.09585209 0.59289670

#######################
#First Rucker Analysis#
#######################

mr_rucker(dat_1_pre)

#[[1]]$intercept
#Method    Estimate          SE       CI_low      CI_upp          P
#1  Egger fixed effects 0.004958444 0.001931704 0.0011723735 0.008744515 0.01026195
#2 Egger random effects 0.004958444 0.002424228 0.0002070455 0.009709843 0.02040916

#[[1]]$Q
#Method         Q df            P
#1   Q_ivw 105.81041 64 0.0007887152
#2 Q_egger  99.22156 63 0.0024400684
#3  Q_diff   6.58885  1 0.0102619464

#[[1]]$res
#[1] "D"

#[[1]]$selected
#Method nsnp    Estimate         SE    CI_low     CI_upp         P
#4 Rucker   65 -0.09695645 0.08848203 -0.270378 0.07646513 0.2773482

#############
#CONCLUSIONS#
#############

#I FOUND NO EVIDENCE FOR HETEROGENEITY WITH ANY TEST.
#LET'S MAKE THE FIGURES.

dat_1_pre$labels <- NA

tiff("BMI_sedentary_original.png", units="in", width=10, height=10, res=300)
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

saveRDS(dat_1_post, "BMI_ST_MR_DF_after_outlier_extraction")

##########################################
#Calculating standard QC for the variants#
##########################################

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)
mF  = mean(F_)

print(mF)

#57.01487

#Cannot find the code for Isq, but we will find it.

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)

#0.9827584

##########################
#First heterogeneity test#
##########################

library(meta)

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_post$mr, dat_1_post$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0002 [0.0000; 0.0199]; tau = 0.0143 [0.0000; 0.1411];
#I^2 = 0.5% [0.0%; 31.5%]; H = 1.00 [1.00; 1.21]

#Test of heterogeneity:
#  Q d.f. p-value
#56.29   56  0.4640

#####################
#Generating analysis#
#####################

mr(dat_1_post)

#id.exposure id.outcome   outcome exposure                    method nsnp           b         se       pval
#1         BMI  sedentary sedentary      BMI                  MR Egger   57 -0.02371395 0.07583585 0.75569037
#2         BMI  sedentary sedentary      BMI           Weighted median   57  0.02831167 0.04020368 0.48130445
#3         BMI  sedentary sedentary      BMI Inverse variance weighted   57  0.06524322 0.02652256 0.01389702
#4         BMI  sedentary sedentary      BMI               Simple mode   57 -0.06272426 0.10343729 0.54670032
#5         BMI  sedentary sedentary      BMI             Weighted mode   57 -0.06272426 0.08369548 0.45673228

####################
#CALCULATING RUCKER#
####################

mr_rucker(dat_1_post)

#[[1]]$intercept
#Method    Estimate          SE       CI_low
#1  Egger fixed effects 0.002596212 0.002039444 -0.001401025
#2 Egger random effects 0.002596212 0.002074938 -0.001470593
#CI_upp         P
#1 0.006593449 0.2030182
#2 0.006663016 0.1054265

#[[1]]$Q
#Method         Q df         P
#1   Q_ivw 58.551614 56 0.3819552
#2 Q_egger 56.931086 55 0.4030557
#3  Q_diff  1.620528  1 0.2030182

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp   Estimate         SE     CI_low    CI_upp          P
#1 Rucker   57 0.06524322 0.02593821 0.01440526 0.1160812 0.01189202

#We are gonna double-check:

mr_rucker_jackknife(dat_1_post)

#We checked with a jackknife and it seems that 
#the model for IVW is the best.

#################
#PRINTING IMAGES#
#################

dat_1_post$labels <- NA

tiff("BMI_sedentary_post.png", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()