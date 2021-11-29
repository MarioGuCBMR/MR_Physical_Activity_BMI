##############
#INTRODUCTION#
##############

#This is a code to run 2SMR with sedentary as exposure and BMI from Locke et al as outcome.

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

#######################
#Loading original data#
#######################

memory.limit(size=80000000)

sedentary <- fread("~/CURATED_DATA/Sedentary_doherty.txt")

BMI <- fread("~/CURATED_DATA/BMI_Locke_Chr_Pos.txt")

####################################################
#Getting genome-wide significant SNPs for sedentary#
####################################################

sedentary_gw <- sedentary[which(as.numeric(sedentary$P_BOLT_LMM_INF) < 0.00000005),]

summary(sedentary_gw$INFO) #good.
summary(sedentary_gw$A1FREQ) #good.
summary(sedentary_gw$P_BOLT_LMM_INF) #good.

#Getting them independent:

sedentary_gw_ind <- sedentary_gw

sedentary_gw_ind$rsid <- sedentary_gw_ind$SNP
sedentary_gw_ind$pval <- sedentary_gw_ind$P_BOLT_LMM_INF

sedentary_gw_ind <- ieugwasr::ld_clump(sedentary_gw_ind) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with just two SNPs..., which is less than the minimum for most 2SMR analysis.
#We are gonna relax the filter.

#################################################################
#Getting P < 5E-07 as genome-wide significant SNPs for sedentary#
#################################################################

sedentary_gw_07 <- sedentary[which(as.numeric(sedentary$P_BOLT_LMM_INF) < 0.0000005),]

summary(sedentary_gw_07$INFO) #good.
summary(sedentary_gw_07$A1FREQ) #good.
summary(sedentary_gw_07$P_BOLT_LMM_INF) #good.

#We make sure to get those legal alleles:

yes_vect <- c("A", "G", "T", "C")

sedentary_gw_07 <- sedentary_gw_07[which(sedentary_gw_07$ALLELE1%in%yes_vect),]
sedentary_gw_07 <- sedentary_gw_07[which(sedentary_gw_07$ALLELE0%in%yes_vect),]

#Let's get the lead, independent variants:

sedentary_gw_ind_07 <- sedentary_gw_07

which(sedentary_gw_ind_07$CHR == 6) #no chromosome 6. No need to remove MHC variants.
 
sedentary_gw_ind_07$rsid <- sedentary_gw_ind_07$SNP
sedentary_gw_ind_07$pval <- sedentary_gw_ind_07$P_BOLT_LMM_INF

sedentary_gw_ind_07 <- ieugwasr::ld_clump_local(sedentary_gw_ind_07, bfile = "~/1K/1kg.v3/1kg.v3/EUR", plink_bin ="~/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

###########################
#Matching with the outcome#
###########################

#We are gonna match with the GIANT RSID, the chromosome and position and with the newest RSID assigned to that position (if exists.)

BMI_oldrsid_2_sedentary <- BMI[which(BMI$SNP%in%sedentary_gw_ind_07$SNP),] #6/14 as the classic 2SMR versions.
BMI_newrsid_2_sedentary <- BMI[which(BMI$newRSID%in%sedentary_gw_ind_07$SNP),] #6/14 same as above.
BMI_chrpos_2_sedentary <- BMI[which(BMI$chr_pos%in%sedentary_gw_ind_07$chr_pos),] #6/14 same as above.

#Good!! Any of these work just fine.
#We are going to set the variants for later. 
#The pre proxy dataframe is really important because in the clumping we are taking into account the variants
#that are only matching with the BMI. 

sedentary_BMI_pre_proxy <- sedentary_gw_ind_07[which(sedentary_gw_ind_07$SNP%in%BMI_chrpos_2_sedentary$SNP),]
BMI_sedentary_pre_proxy <- BMI_chrpos_2_sedentary

#################################################################
#Getting proxies for the SNPs that we do not have in the outcome#
#################################################################

missings_snps <- sedentary_gw_ind_07$SNP[which(!(sedentary_gw_ind_07$SNP%in%BMI_chrpos_2_sedentary$SNP))]

#Get the proxies:

proxies <- ld_proxies(missings_snps) #do not worry about the error. It does not affect obtaining the proxies.

proxies_clean <- as.character(unlist(proxies$variation2)) #196 proxies.

#Since we have only RSIDs we are gonna make a trick and check if they are
#the same in BMI for old or new RSIDs. 
#We are working with BMI first, because it has less variants and it becomes a bottleneck.

BMI_oldrsid_2_proxies <- BMI[which(BMI$SNP%in%proxies_clean),] #91
BMI_newrsid_2_proxies <- BMI[which(BMI$newRSID%in%proxies_clean),] #91

#We obtain 91 in each case. Are they the same with old and new RSIDs?

which(BMI_oldrsid_2_proxies$SNP != BMI_oldrsid_2_proxies$newRSID) #yup.

#And now we match the sedentary with this:

sedentary_BMI_post_proxy <- sedentary[which(sedentary$SNP%in%BMI_oldrsid_2_proxies$SNP),]
BMI_sedentary_post_proxy <- BMI_oldrsid_2_proxies[which(BMI_oldrsid_2_proxies$SNP%in%sedentary_BMI_post_proxy$SNP),]

############################
#Merging pre and post proxy#
############################

sedentary_BMI_post_proxy$rsid <- sedentary_BMI_post_proxy$SNP
sedentary_BMI_post_proxy$pval <- sedentary_BMI_post_proxy$P_BOLT_LMM_INF
sedentary_BMI_post_proxy$id <- "sedentary"
sedentary_BMI_pre_proxy$id <- "sedentary"

sedentary_pre_merge <- rbind(sedentary_BMI_pre_proxy, sedentary_BMI_post_proxy)

#Now the same for BMI:

BMI_pre_merge <- rbind(BMI_sedentary_pre_proxy, BMI_sedentary_post_proxy)

BMI_pre_merge$id <- "BMI"

###########################################################
#Before merging we need to do the pruning for the exposure#
###########################################################

sedentary_final_set <- ieugwasr::ld_clump_local(sedentary_pre_merge, bfile = "~/1K/1kg.v3/1kg.v3/EUR", plink_bin ="~/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with 11 SNPs!

################################
#Preparing data for TwoSampleMR#
################################

BMI_final_set <- BMI_pre_merge[which(BMI_pre_merge$SNP%in%sedentary_final_set$SNP),]

#We should have the same SNPs in each.
#Perfect, now we get the columns for TwoSampleMR

sedentary_4_merge <- sedentary_final_set %>%
  select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF, id)

colnames(sedentary_4_merge) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "id.exposure")

BMI_4_merge <- BMI_final_set %>%
  select(SNP, A1, A2, Freq1.Hapmap, b, se, p, id, N)

colnames(BMI_4_merge) <- c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "id.outcome", "samplesize.outcome")

#######################
#Finally we merge them#
#######################

sedentary_4_merge$exposure <- "ST"
BMI_4_merge$outcome <- "BMI"

dat_1_pre <- harmonise_data(sedentary_4_merge, BMI_4_merge, action = 3)

dat_1_pre <- dat_1_pre[which(dat_1_pre$palindromic == FALSE),]
dat_1_pre <- dat_1_pre[which(dat_1_pre$remove == FALSE),]

######################
#Let's check the data#
######################

summary(dat_1_pre$eaf.exposure) #perfect
summary(dat_1_pre$eaf.outcome) #we need to remove one.

dat_1_pre <- dat_1_pre[which(is.na(dat_1_pre$eaf.outcome) == FALSE),] #we end up with 8 SNPs.

###############################
#Now let's add the sample size#
###############################

dat_1_pre$samplesize.exposure <- 91105
dat_1_pre$exposure <- "sedentary"
dat_1_pre$outcome <- "BMI"

#Now we can perform steiger:

dat_1_pre <- steiger_filtering(dat_1_pre)

dat_1_pre <- dat_1_pre[which(dat_1_pre$steiger_dir == TRUE),]

saveRDS(dat_1_pre, "ST_BMI_MR_DF")

##############################
#THUS we can run our analysis#
##############################

###############################################
#Checking mF and Isq and chromose and position#
###############################################

sedentary_pos_test <- sedentary[which(sedentary$SNP%in%dat_1_pre$SNP),]

which(sedentary_pos_test$CHR == 6) #no worries with this one.

F_ = ((dat_1_pre$beta.exposure)^2)/((dat_1_pre$se.exposure)^2)
mF  = mean(F_)

print(mF)
# 29.37228

Isq(dat_1_pre$beta.exposure, dat_1_pre$se.exposure)
#0.9609034

##########################
#First heterogeneity test#
##########################

library(meta)

dat_1_pre$mr <- dat_1_pre$beta.outcome/dat_1_pre$beta.exposure
dat_1_pre$mr_se <- ((dat_1_pre$mr*((dat_1_pre$se.exposure/dat_1_pre$beta.exposure)^2+(dat_1_pre$se.outcome/dat_1_pre$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_pre$mr, dat_1_pre$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0589 [0.0112; 0.3300]; tau = 0.2427 [0.1060; 0.5745];
#I^2 = 69.7% [36.8%; 85.4%]; H = 1.82 [1.26; 2.62]

#Test of heterogeneity:
#  Q d.f. p-value
#23.08    7  0.0017

###################
#First MR analysis#
###################

mr(dat_1_pre)

#id.exposure id.outcome outcome  exposure                    method nsnp           b         se      pval
#1   sedentary        BMI     BMI sedentary                  MR Egger    8  0.46163922 0.64266166 0.4995608
#2   sedentary        BMI     BMI sedentary           Weighted median    8  0.02119388 0.09444378 0.8224403
#3   sedentary        BMI     BMI sedentary Inverse variance weighted    8  0.03836607 0.10760804 0.7214397
#4   sedentary        BMI     BMI sedentary               Simple mode    8  0.11506251 0.18651897 0.5568284
#5   sedentary        BMI     BMI sedentary             Weighted mode    8 -0.05418201 0.18456499 0.7776024

#######################
#First Rucker Analysis#
#######################

mr_rucker(dat_1_pre)

#[[1]]$intercept
#Method    Estimate         SE      CI_low      CI_upp
#1  Egger fixed effects -0.01204261 0.00859454 -0.02888760 0.004802379
#2 Egger random effects -0.01204261 0.01800403 -0.04732987 0.023244644
#P
#1 0.1611563
#2 0.7482153

#[[1]]$Q
#Method         Q df            P
#1   Q_ivw 28.293017  7 0.0001945649
#2 Q_egger 26.329674  6 0.0001932574
#3  Q_diff  1.963342  1 0.1611563321

#[[1]]$res
#[1] "B"

#[[1]]$selected
#Method nsnp   Estimate       SE     CI_low   CI_upp         P
#2 Rucker    8 0.03836607 0.107608 -0.1725418 0.249274 0.7214397

#############
#CONCLUSIONS#
#############

dat_1_pre$labels <- NA

tiff("Sedentary_BMI_Proxy_original_P07.png", units="in", width=10, height=10, res=300)
mr_plots(dat_1_pre)
dev.off()

##########################
#Let's remove the outlier#
##########################

##################
#REMOVING OUTLIER#
##################

library(RadialMR)

radial_input <- format_radial(dat_1_pre$beta.exposure, dat_1_pre$beta.outcome, dat_1_pre$se.exposure, dat_1_pre$se.outcome, dat_1_pre$SNP)

radial_output <- ivw_radial(radial_input, alpha = 0.05, weights = 3, tol = 0.0001) #mF from RadialMR is 0.13/1 good.

outliers <- radial_output$outliers$SNP 

dat_1_post <- dat_1_pre[-which(dat_1_pre$SNP%in%outliers),] 

saveRDS(dat_1_post, "ST_BMI_MR_DF_after_outlier_extraction")

#######################
#Analysis post outlier#
#######################

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)
mF  = mean(F_)

print(mF)
#27.97988

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#0.9607653

#And now we check out once again:

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_post$mr, dat_1_post$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0046 [0.0000; 0.2269]; tau = 0.0682 [0.0000; 0.4763];
#I^2 = 15.0% [0.0%; 82.3%]; H = 1.08 [1.00; 2.38]

#Test of heterogeneity:
#  Q d.f. p-value
#4.70    4  0.3191

##################
#Quantifying data#
##################

mr(dat_1_post)

#id.exposure id.outcome outcome  exposure                    method nsnp          b         se      pval
#1   sedentary        BMI     BMI sedentary                  MR Egger    5 0.48920463 0.35757764 0.2647400
#2   sedentary        BMI     BMI sedentary           Weighted median    5 0.10132701 0.09544084 0.2883840
#3   sedentary        BMI     BMI sedentary Inverse variance weighted    5 0.11231999 0.07911264 0.1556812
#4   sedentary        BMI     BMI sedentary               Simple mode    5 0.19959321 0.15224172 0.2600459
#5   sedentary        BMI     BMI sedentary             Weighted mode    5 0.08484545 0.14356941 0.5863170

mr_rucker(dat_1_post)

#[[1]]$intercept
#Method    Estimate          SE      CI_low      CI_upp
#1  Egger fixed effects -0.01051476 0.008870840 -0.02790128 0.006871772
#2 Egger random effects -0.01051476 0.009738833 -0.02960252 0.008573006
#P
#1 0.2358922
#2 0.8598561

#[[1]]$Q
#Method        Q df         P
#1   Q_ivw 5.020785  4 0.2851714
#2 Q_egger 3.615810  3 0.3060497
#3  Q_diff 1.404976  1 0.2358922

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp Estimate         SE      CI_low    CI_upp         P
#1 Rucker    5  0.11232 0.07061387 -0.02608067 0.2507206 0.1116947

dat_1_post$labels <- NA

tiff("Sedentary_BMI_Proxy_Post_P07.png", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()
