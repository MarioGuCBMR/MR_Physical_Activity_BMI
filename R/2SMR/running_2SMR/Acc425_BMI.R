##############
#INTRODUCTION#
##############

#This is a code to run 2SMR with Acc425 as exposure and BMI from Locke et al as outcome.

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

acc425 <- fread("~/CURATED_DATA/acc425_klimentidis.txt")

BMI <- fread("~/CURATED_DATA/BMI_Locke_Chr_Pos.txt")

##############################################################
#Getting P < 5E-07 as genome-wide significant SNPs for acc425#
##############################################################

acc425_gw_07 <- acc425[which(as.numeric(acc425$P_BOLT_LMM_INF) < 0.0000005),]

which(acc425_gw_07$CHR == 6)
summary(acc425_gw_07$INFO) #good.
summary(acc425_gw_07$A1FREQ) #good.
summary(acc425_gw_07$P_BOLT_LMM_INF) #good.

#Getting the lead, independent variants:

acc425_gw_ind_07 <- acc425_gw_07

acc425_gw_ind_07$rsid <- acc425_gw_ind_07$SNP
acc425_gw_ind_07$pval <- acc425_gw_ind_07$P_BOLT_LMM_INF

acc425_gw_ind_07 <- ieugwasr::ld_clump_local(acc425_gw_ind_07, bfile = "~1K/1kg.v3/1kg.v3/EUR", plink_bin ="~/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

###########################
#Matching with the outcome#
###########################

#We are going to do the matching with several sources of info:

BMI_oldrsid_2_acc425 <- BMI[which(BMI$SNP%in%acc425_gw_ind_07$SNP),] #2/8 as the classic 2SMR versions.
BMI_newrsid_2_acc425 <- BMI[which(BMI$newRSID%in%acc425_gw_ind_07$SNP),] #2/8 same as above.
BMI_chrpos_2_acc425 <- BMI[which(BMI$chr_pos%in%acc425_gw_ind_07$chr_pos),] #2/8 same as above.

#Good!! Any of these work just fine.

Acc425_BMI_pre_proxy <- acc425_gw_ind_07[which(acc425_gw_ind_07$SNP%in%BMI_chrpos_2_acc425$SNP),]
BMI_acc425_pre_proxy <- BMI_chrpos_2_acc425

#################################################################
#Getting proxies for the SNPs that we do not have in the outcome#
#################################################################

missings_snps <- acc425_gw_ind_07$SNP[which(!(acc425_gw_ind_07$SNP%in%BMI_chrpos_2_acc425$SNP))]

#Get the proxies:

proxies <- ld_proxies(missings_snps)

proxies_clean <- as.character(unlist(proxies$variation2)) #2254 proxies.

#We are going to match the proxies in BMI first since it has less variants and 
#will works as a bottleneck.
#We will need to check the new and the old RSIDS.

BMI_oldrsid_2_proxies <- BMI[which(BMI$SNP%in%proxies_clean),] #324
BMI_newrsid_2_proxies <- BMI[which(BMI$newRSID%in%proxies_clean),] #324

#Do we have the same variants?

which(BMI_oldrsid_2_proxies$SNP != BMI_oldrsid_2_proxies$newRSID) #yup.

#And now we match the Acc425 with this:

acc425_BMI_post_proxy <- acc425[which(acc425$SNP%in%BMI_oldrsid_2_proxies$SNP),]
BMI_acc425_post_proxy <- BMI_oldrsid_2_proxies[which(BMI_oldrsid_2_proxies$SNP%in%acc425_BMI_post_proxy$SNP),]

################################################################################
#We got the proxies in the bottleneck way, let's check, just in case, the other#
################################################################################

acc425_proxy_trial <- acc425[which(acc425$SNP%in%proxies_clean),] #1242. The best way would be taking those that have a chromosome and position...

#But we have seen that the bottleneck version allows us to see that even in this case it is fine.

BMI_proxy_trial <- BMI[which(BMI$SNP%in%acc425_proxy_trial$SNP),] #the same amount of SNPs.

#The bottleneck works perfectly!!!

#We just need to merge the dataframes from pre and post proxy for each trait.
#Convert them to the TwoSampleMR dataframe
#And run the analysis.

############################
#Merging pre and post proxy#
############################

acc425_BMI_post_proxy$rsid <- acc425_BMI_post_proxy$SNP
acc425_BMI_post_proxy$pval <- acc425_BMI_post_proxy$P_BOLT_LMM_INF
acc425_BMI_post_proxy$id <- "Acc425"
Acc425_BMI_pre_proxy$id <- "Acc425"

acc425_pre_merge <- rbind(Acc425_BMI_pre_proxy, acc425_BMI_post_proxy)

#Now the same for BMI:

BMI_pre_merge <- rbind(BMI_acc425_pre_proxy, BMI_acc425_post_proxy)

BMI_pre_merge$id <- "BMI"

###########################################################
#Before merging we need to do the pruning for the exposure#
###########################################################

acc425_final_set <- ieugwasr::ld_clump_local(acc425_pre_merge, bfile = "~1K/1kg.v3/1kg.v3/EUR", plink_bin ="~/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with 6 SNPs!

################################
#Preparing data for TwoSampleMR#
################################

BMI_final_set <- BMI_pre_merge[which(BMI_pre_merge$SNP%in%acc425_final_set$SNP),]

#We should have the same SNPs in each.
#Perfect, now we get the columns for TwoSampleMR

acc425_4_merge <- acc425_final_set %>%
  select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF, id)

colnames(acc425_4_merge) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "id.exposure")

BMI_4_merge <- BMI_final_set %>%
  select(SNP, A1, A2, Freq1.Hapmap, b, se, p, id, N)

colnames(BMI_4_merge) <- c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "id.outcome", "samplesize.outcome")

#######################
#Finally we merge them#
#######################

acc425_4_merge$exposure <- "Acc425"
BMI_4_merge$outcome <- "BMI"

dat_1_pre <- harmonise_data(acc425_4_merge, BMI_4_merge, action = 3)

dat_1_pre <- dat_1_pre[which(dat_1_pre$palindromic == FALSE),]
dat_1_pre <- dat_1_pre[which(dat_1_pre$remove == FALSE),]

######################
#Let's check the data#
######################

summary(dat_1_pre$eaf.exposure) #perfect
summary(dat_1_pre$eaf.outcome) #perfect.

###############################
#Now let's add the sample size#
###############################

dat_1_pre$samplesize.exposure <- 90667
dat_1_pre$exposure <- "Acc425"
dat_1_pre$outcome <- "BMI"

#Now we can perform steiger:

dat_1_pre <- steiger_filtering(dat_1_pre)

dat_1_pre <- dat_1_pre[which(dat_1_pre$steiger_dir == TRUE),]

##############################
#THUS we can run our analysis#
##############################

saveRDS(dat_1_pre, "Acc425_BMI_MR_DF")

###############################################
#Checking mF and Isq and chromose and position#
###############################################

Acc425_pos_test <- acc425[which(acc425$SNP%in%dat_1_pre$SNP),]

which(Acc425_pos_test$CHR == 6) #no worries with this one.

F_ = ((dat_1_pre$beta.exposure)^2)/((dat_1_pre$se.exposure)^2)
mF  = mean(F_)

print(mF)

#28.65615

Isq(dat_1_pre$beta.exposure, dat_1_pre$se.exposure)

#0.9717361

##########################
#First heterogeneity test#
##########################

library(meta)

dat_1_pre$mr <- dat_1_pre$beta.outcome/dat_1_pre$beta.exposure
dat_1_pre$mr_se <- ((dat_1_pre$mr*((dat_1_pre$se.exposure/dat_1_pre$beta.exposure)^2+(dat_1_pre$se.outcome/dat_1_pre$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_pre$mr, dat_1_pre$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0 [0.0000; 0.1422]; tau = 0 [0.0000; 0.3771];
#I^2 = 0.0% [0.0%; 49.5%]; H = 1.00 [1.00; 1.41]

#Test of heterogeneity:
#  Q d.f. p-value
#1.65    4  0.8003

###################
#First MR analysis#
###################

mr(dat_1_pre)

#id.exposure id.outcome outcome exposure                    method nsnp          b         se       pval
#1      Acc425        BMI     BMI   Acc425                  MR Egger    5  1.3253160 2.20819691 0.59069563
#2      Acc425        BMI     BMI   Acc425           Weighted median    5 -0.1762778 0.10060734 0.07975038
#3      Acc425        BMI     BMI   Acc425 Inverse variance weighted    5 -0.1737542 0.08451525 0.03979292
#4      Acc425        BMI     BMI   Acc425               Simple mode    5 -0.2025959 0.12600868 0.18316020
#5      Acc425        BMI     BMI   Acc425             Weighted mode    5 -0.1871664 0.12030467 0.19474392

#######################
#First Rucker Analysis#
#######################

mr_rucker(dat_1_pre)

#[[1]]$intercept
#Method    Estimate         SE     CI_low     CI_upp
#1  Egger fixed effects -0.03841132 0.03663385 -0.1102124 0.03338971
#2 Egger random effects -0.03841132 0.03663385 -0.1102124 0.03338971
#P
#1 0.2943992
#2 0.8528004

#[[1]]$Q
#Method         Q df         P
#1   Q_ivw 1.7209625  4 0.7869063
#2 Q_egger 1.2594271  3 0.7387893
#3  Q_diff 0.4615354  1 0.4969073

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp   Estimate         SE     CI_low      CI_upp           P
#1 Rucker    5 -0.1737542 0.05543586 -0.2824065 -0.06510195 0.001722475

#############
#CONCLUSIONS#
#############

#I FOUND NO EVIDENCE FOR HETEROGENEITY WITH ANY TEST.
#LET'S MAKE THE FIGURES.

dat_1_pre$labels <- NA

tiff("Acc455_BMI_Proxy_original.png", units="in", width=10, height=10, res=300)
mr_plots(dat_1_pre)
dev.off()
