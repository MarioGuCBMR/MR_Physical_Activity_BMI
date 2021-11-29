##############
#INTRODUCTION#
##############

#This is a code to run 2SMR with BMI as exposure and moderate physical activity as outcome.

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

###################
#Loading functions#
###################

mr_plots <- function(dat)
{
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

moderate <- fread("~/CURATED_DATA/moderate_doherty.txt")

BMI <- fread("~/CURATED_DATA/BMI_Locke_Chr_Pos.txt")

###################################################
#Getting genome-wide significant SNPs for moderate#
###################################################

BMI_gw <- BMI[which(as.numeric(BMI$p) < 0.00000005),]

summary(BMI_gw$Freq1.Hapmap) #good. #the NAs will be removed.
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

#All good.

#####################################
#Checking if we have merged variants#
#####################################

#The GIANT data is a bit old, so some rsIDs might be outdated and won't be found in the 1000G reference panel.
#For this reason, in the curation of the data we included new rsIDs for the same chr_pos. 
#Now I am going to check whether there are any discrepancies between GIANT rsIDS and those found in the merged rsIDs database.

which(BMI_gw_ind_oldRSID$SNP != BMI_gw_ind_oldRSID$newRSID) #we have one different!

#We will have to check if we have different results with the old and the new rsIDs.

################################
#LD CLUMPING WITH ORIGINAL RSID#
################################

BMI_gw_ind_oldRSID_ <- ieugwasr::ld_clump_local(BMI_gw_ind_oldRSID, bfile = "~/1K/1kg.v3/1kg.v3/EUR", plink_bin ="~/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We got 68 SNPs, let's check with the new RSIDs

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

moderate_oldrsid_2_BMI <- moderate[which(moderate$chr_pos%in%BMI_gw_ind_oldRSID_$chr_pos),] #68/68 as the classic 2SMR versions.
moderate_newrsid_2_BMI <- moderate[which(moderate$chr_pos%in%BMI_gw_ind_newRSID_$chr_pos),] #68/68 but one SNP has changed.

#The lead SNP changes taking into account which version we are using. 
#We will have to compare, though of course, I doubt the results change a lot since they are bound to be in the same loci.

################################
#Preparing data for TwoSampleMR#
################################

#We are going to do this twice to check whether we wind up with any differences or not:

moderate_old_set <- moderate_oldrsid_2_BMI
BMI_old_set <- BMI_gw_ind_oldRSID_

moderate_old_set$id <- "moderate"
BMI_old_set$id <- "BMI"

##############################
#FIRST WITH THE ORIGINAL DATA#
##############################

#We should have the same SNPs in each.
#Perfect, now we get the columns for TwoSampleMR

moderate_4_merge_old <- moderate_old_set %>%
  select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF, id)

colnames(moderate_4_merge_old) <- c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "id.outcome")

BMI_4_merge_old <- BMI_old_set %>%
  select(SNP, A1, A2, Freq1.Hapmap, b, se, p, id, N)

colnames(BMI_4_merge_old) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "id.exposure", "samplesize.exposure")

#######################
#Finally we merge them#
#######################

BMI_4_merge_old$exposure <- "BMI"
moderate_4_merge_old$outcome <- "Moderate"

dat_1_pre_old <- harmonise_data(BMI_4_merge_old, moderate_4_merge_old, action = 3)

dat_1_pre_old <- dat_1_pre_old[which(dat_1_pre_old$palindromic == FALSE),]
dat_1_pre_old <- dat_1_pre_old[which(dat_1_pre_old$remove == FALSE),]

#########################################
#WE REPEAT THIS BUT WITH NEWRSID VERSION#
#########################################

#We are going to do this twice to check whether we wind up with any differences or not:

moderate_new_set <- moderate_newrsid_2_BMI
BMI_new_set <- BMI_gw_ind_newRSID_

moderate_new_set$id <- "moderate"
BMI_new_set$id <- "BMI"

#######################
#NOW WITH THE NEW DATA#
#######################

#We should have the same SNPs in each.
#Perfect, now we get the columns for TwoSampleMR

moderate_4_merge_new <- moderate_new_set %>%
  select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF, id)

colnames(moderate_4_merge_new) <- c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "id.outcome")

BMI_4_merge_new <- BMI_new_set %>%
  select(newRSID, A1, A2, Freq1.Hapmap, b, se, p, id, N)

colnames(BMI_4_merge_new) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "id.exposure", "samplesize.exposure")

#######################
#Finally we merge them#
#######################

BMI_4_merge_new$exposure <- "BMI"
moderate_4_merge_new$outcome <- "Moderate"

dat_1_pre_new <- harmonise_data(BMI_4_merge_new, moderate_4_merge_new, action = 3)

dat_1_pre_new <- dat_1_pre_new[which(dat_1_pre_new$palindromic == FALSE),]
dat_1_pre_new <- dat_1_pre_new[which(dat_1_pre_new$remove == FALSE),]

#We get one more SNP than before.
#Hence, we stick with the new one:

dat_1_pre <- dat_1_pre_new

######################
#Let's check the data#
######################

summary(dat_1_pre$eaf.exposure) #perfect
summary(dat_1_pre$eaf.outcome) #perfect.

###############################
#Now let's add the sample size#
###############################

dat_1_pre$samplesize.outcome <- 91105
dat_1_pre$outcome <- "moderate"
dat_1_pre$exposure <- "BMI"

#Now we can perform steiger:

dat_1_pre <- steiger_filtering(dat_1_pre)

dat_1_pre <- dat_1_pre[which(dat_1_pre$steiger_dir == TRUE),] 

##############################
#THUS we can run our analysis#
##############################

saveRDS(dat_1_pre, "BMI_Moderate_MR_DF")

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
#  tau^2 = 0.0121 [0.0024; 0.0519]; tau = 0.1099 [0.0494; 0.2278];
#I^2 = 23.6% [0.0%; 44.2%]; H = 1.14 [1.00; 1.34]

#Test of heterogeneity:
#  Q d.f. p-value
#83.73   64  0.0496

###################
#First MR analysis#
###################

mr(dat_1_pre)

#id.exposure id.outcome  outcome exposure                    method nsnp            b         se       pval
#1         BMI   moderate moderate      BMI                  MR Egger   65  0.170135830 0.08081788 0.03926953
#2         BMI   moderate moderate      BMI           Weighted median   65  0.003900416 0.03752978 0.91722603
#3         BMI   moderate moderate      BMI Inverse variance weighted   65 -0.049119716 0.02914157 0.09188160
#4         BMI   moderate moderate      BMI               Simple mode   65 -0.003585929 0.06676111 0.95733126
#5         BMI   moderate moderate      BMI             Weighted mode   65  0.017994983 0.05043587 0.72242319

#######################
#First Rucker Analysis#
#######################

mr_rucker(dat_1_pre)

#[[1]]$intercept
#Method     Estimate          SE      CI_low
#1  Egger fixed effects -0.006391344 0.001939194 -0.01019209
#2 Egger random effects -0.006391344 0.002214186 -0.01073107
#CI_upp            P
#1 -0.002590594 0.0009811492
#2 -0.002051619 0.9980525056

#[[1]]$Q
#Method        Q df            P
#1   Q_ivw 92.99745 64 0.0104109761
#2 Q_egger 82.13464 63 0.0531066124
#3  Q_diff 10.86281  1 0.0009811492

#[[1]]$res
#[1] "C"

#[[1]]$selected
#Method nsnp  Estimate         SE    CI_low    CI_upp          P
#3 Rucker   65 0.1701358 0.07078065 0.0314083 0.3088634 0.01918428

#############
#CONCLUSIONS#
#############

#I FOUND NO EVIDENCE FOR HETEROGENEITY WITH ANY TEST.
#LET'S MAKE THE FIGURES.

dat_1_pre$labels <- NA

tiff("BMI_moderate_original.png", units="in", width=10, height=10, res=300)
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

saveRDS(dat_1_post, "BMI_Moderate_MR_DF_after_outlier_extraction")

##########################################
#Calculating standard QC for the variants#
##########################################

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)
mF  = mean(F_)

print(mF)

#61.79095

#Cannot find the code for Isq, but we will find it.

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)

#0.9840239

##########################
#First heterogeneity test#
##########################

library(meta)

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_post$mr, dat_1_post$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0 [0.0000; 0.0084]; tau = 0 [0.0000; 0.0918];
#I^2 = 0.0% [0.0%; 8.9%]; H = 1.00 [1.00; 1.05]

#Test of heterogeneity:
#  Q d.f. p-value
#40.51   54  0.9131

#####################
#Generating analysis#
#####################

mr(dat_1_post)

#id.exposure id.outcome  outcome exposure                    method nsnp            b         se       pval
#1         BMI   moderate moderate      BMI                  MR Egger   55  0.159638962 0.07410200 0.03578081
#2         BMI   moderate moderate      BMI           Weighted median   55  0.002951477 0.03611421 0.93486444
#3         BMI   moderate moderate      BMI Inverse variance weighted   55 -0.055927261 0.02554654 0.02858010
#4         BMI   moderate moderate      BMI               Simple mode   55  0.007233163 0.06889888 0.91677885
#5         BMI   moderate moderate      BMI             Weighted mode   55  0.024308112 0.05786817 0.67610828

####################
#CALCULATING RUCKER#
####################

mr_rucker(dat_1_post)

#[[1]]$intercept
#Method     Estimate          SE       CI_low
#1  Egger fixed effects -0.006439501 0.001646064 -0.009665727
#2 Egger random effects -0.006439501 0.001646064 -0.009665727
#CI_upp            P
#1 -0.003213275 9.151222e-05
#2 -0.003213275 9.999542e-01

###############################
#WE RUN RADIALMR ONE MORE TIME#
###############################

#[[1]]$Q
#Method         Q df           P
#1   Q_ivw 42.863622 54 0.862199237
#2 Q_egger 33.259617 53 0.984569526
#3  Q_diff  9.604005  1 0.001941535

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp    Estimate         SE     CI_low      CI_upp          P
#1 Rucker   55 -0.05592726 0.02276039 -0.1005368 -0.01131772 0.01400172

#We are gonna double-check:

mr_rucker_jackknife(dat_1_post)

#We checked with a jackknife and it seems that 
#the model for IVW is the best.

#################
#PRINTING IMAGES#
#################

dat_1_post$labels <- NA

tiff("BMI_moderate_post.png", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()

###############################
#WE RUN RADIALMR ONE MORE TIME#
###############################

library(RadialMR)

radial_input <- format_radial(dat_1_post$beta.exposure, dat_1_post$beta.outcome, dat_1_post$se.exposure, dat_1_post$se.outcome, dat_1_post$SNP)

radial_output <- ivw_radial(radial_input, 0.05, 3)  #no outliers.

outliers <- radial_output$outliers$SNP 

dat_1_post <- dat_1_pre[which(!(dat_1_pre$SNP%in%outliers)),]

#NO MORE OUTLIERS

