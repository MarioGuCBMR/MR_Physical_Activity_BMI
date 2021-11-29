##############
#INTRODUCTION#
##############

#This is a code to run 2SMR with moderate as exposure and BMI from Locke et al as outcome.

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

Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

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

moderate <- fread("~/CURATED_DATA/Moderate_doherty.txt")

BMI <- fread("~/CURATED_DATA/BMI_Locke_Chr_Pos.txt")

###################################################
#Getting genome-wide significant SNPs for moderate#
###################################################

#################################################################
#Getting P < 5E-07 as genome-wide significant SNPs for moderate#
#################################################################

moderate_gw_07 <- moderate[which(as.numeric(moderate$P_BOLT_LMM_INF) < 0.0000005),]

summary(moderate_gw_07$INFO) #good.
summary(moderate_gw_07$A1FREQ) #good.
summary(moderate_gw_07$P_BOLT_LMM_INF) #good.

#We make sure to get those legal alleles:

yes_vect <- c("A", "G", "T", "C")

moderate_gw_07 <- moderate_gw_07[which(moderate_gw_07$ALLELE1%in%yes_vect),]
moderate_gw_07 <- moderate_gw_07[which(moderate_gw_07$ALLELE0%in%yes_vect),]

moderate_gw_07 <- moderate_gw_07[order(moderate_gw_07$SNP),]

head(moderate_gw_07)
tail(moderate_gw_07)

#Getting them independent:

moderate_gw_ind_07 <- moderate_gw_07

which(moderate_gw_ind_07$CHR == 6) #no chromosome 6.

moderate_gw_ind_07$rsid <- moderate_gw_ind_07$SNP
moderate_gw_ind_07$pval <- moderate_gw_ind_07$P_BOLT_LMM_INF

moderate_gw_ind_07 <- ieugwasr::ld_clump_local(moderate_gw_ind_07, bfile = "~/1K/1kg.v3/1kg.v3/EUR", plink_bin ="~/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

###########################
#Matching with the outcome#
###########################

#We are going to run it with chromosome and position (build 37), the original GIANT RSID and the newest available
#RSIDs, if there are any, for the chromsome and position.

BMI_oldrsid_2_moderate <- BMI[which(BMI$SNP%in%moderate_gw_ind_07$SNP),] #1/4 as the classic 2SMR versions.
BMI_newrsid_2_moderate <- BMI[which(BMI$newRSID%in%moderate_gw_ind_07$SNP),] #1/4 same as above.
BMI_chrpos_2_moderate <- BMI[which(BMI$chr_pos%in%moderate_gw_ind_07$chr_pos),] #1/4 same as above.

#Good!! Any of these work just fine.

moderate_BMI_pre_proxy <- moderate_gw_ind_07[which(moderate_gw_ind_07$SNP%in%BMI_chrpos_2_moderate$SNP),]
BMI_moderate_pre_proxy <- BMI_chrpos_2_moderate

#################################################################
#Getting proxies for the SNPs that we do not have in the outcome#
#################################################################

missings_snps <- moderate_gw_ind_07$SNP[which(!(moderate_gw_ind_07$SNP%in%BMI_chrpos_2_moderate$SNP))]

#Get the proxies:

proxies <- ld_proxies(missings_snps)

proxies_clean <- as.character(unlist(proxies$variation2)) #110 proxies.

#We are going to match the proxies in BMI first since it has less variants and 
#will works as a bottleneck.

BMI_oldrsid_2_proxies <- BMI[which(BMI$SNP%in%proxies_clean),] #40
BMI_newrsid_2_proxies <- BMI[which(BMI$newRSID%in%proxies_clean),] #40

#Are the GIANT SNPs the latest available for the chromosome and position in the proxies?

which(BMI_oldrsid_2_proxies$SNP != BMI_oldrsid_2_proxies$newRSID) #yup.

#And now we match the moderate with this:

moderate_BMI_post_proxy <- moderate[which(moderate$SNP%in%BMI_oldrsid_2_proxies$SNP),]
BMI_moderate_post_proxy <- BMI_oldrsid_2_proxies[which(BMI_oldrsid_2_proxies$SNP%in%moderate_BMI_post_proxy$SNP),]

############################
#Merging pre and post proxy#
############################

moderate_BMI_post_proxy$rsid <- moderate_BMI_post_proxy$SNP
moderate_BMI_post_proxy$pval <- moderate_BMI_post_proxy$P_BOLT_LMM_INF
moderate_BMI_post_proxy$id <- "moderate"
moderate_BMI_pre_proxy$id <- "moderate"

moderate_pre_merge <- rbind(moderate_BMI_pre_proxy, moderate_BMI_post_proxy)

#Now the same for BMI:

BMI_pre_merge <- rbind(BMI_moderate_pre_proxy, BMI_moderate_post_proxy)

BMI_pre_merge$id <- "BMI"

###########################################################
#Before merging we need to do the pruning for the exposure#
###########################################################

moderate_final_set <- ieugwasr::ld_clump_local(moderate_pre_merge, bfile = "~/1K/1kg.v3/1kg.v3/EUR", plink_bin ="~/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99) #no need to do anything. LD r2 < 0.001 and kb = 10000 and P = 1,

#We end up with 3 SNPs!

################################
#Preparing data for TwoSampleMR#
################################

BMI_final_set <- BMI_pre_merge[which(BMI_pre_merge$SNP%in%moderate_final_set$SNP),]

#We should have the same SNPs in each.
#Perfect, now we get the columns for TwoSampleMR

moderate_4_merge <- moderate_final_set %>%
  select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P_BOLT_LMM_INF, id)

colnames(moderate_4_merge) <- c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "id.exposure")

BMI_4_merge <- BMI_final_set %>%
  select(SNP, A1, A2, Freq1.Hapmap, b, se, p, id, N)

colnames(BMI_4_merge) <- c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "id.outcome", "samplesize.outcome")

#######################
#Finally we merge them#
#######################

moderate_4_merge$exposure <- "Moderate"
BMI_4_merge$outcome <- "BMI"

dat_1_pre <- harmonise_data(moderate_4_merge, BMI_4_merge, action = 3)

dat_1_pre <- dat_1_pre[which(dat_1_pre$palindromic == FALSE),]
dat_1_pre <- dat_1_pre[which(dat_1_pre$remove == FALSE),]

######################
#Let's check the data#
######################

summary(dat_1_pre$eaf.exposure) #perfect
summary(dat_1_pre$eaf.outcome) #perfect

dat_1_pre <- dat_1_pre[which(is.na(dat_1_pre$eaf.outcome) == FALSE),] #we end up with 3 SNPs.

###############################
#Now let's add the sample size#
###############################

dat_1_pre$samplesize.exposure <- 91105
dat_1_pre$exposure <- "moderate"
dat_1_pre$outcome <- "BMI"

#Now we can perform steiger:

dat_1_pre <- steiger_filtering(dat_1_pre)

dat_1_pre <- dat_1_pre[which(dat_1_pre$steiger_dir == TRUE),]

##############################
#THUS we can run our analysis#
##############################

saveRDS(dat_1_pre, "Moderate_BMI_MR_DF")

###############################################
#Checking mF and Isq and chromose and position#
###############################################

moderate_pos_test <- moderate[which(moderate$SNP%in%dat_1_pre$SNP),]

which(moderate_pos_test$CHR == 6) #no worries with this one.

F_ = ((dat_1_pre$beta.exposure)^2)/((dat_1_pre$se.exposure)^2)
mF  = mean(F_)

print(mF)
#28.95627

Isq(dat_1_pre$beta.exposure, dat_1_pre$se.exposure)
#0.97673

##########################
#First heterogeneity test#
##########################

library(meta)

dat_1_pre$mr <- dat_1_pre$beta.outcome/dat_1_pre$beta.exposure
dat_1_pre$mr_se <- ((dat_1_pre$mr*((dat_1_pre$se.exposure/dat_1_pre$beta.exposure)^2+(dat_1_pre$se.outcome/dat_1_pre$beta.outcome)^2)^0.5)^2)^0.5
metagen(dat_1_pre$mr, dat_1_pre$mr_se)

#Quantifying heterogeneity:
#  tau^2 = 0.0188 [0.0000; 1.7570]; tau = 0.1371 [0.0000; 1.3255];
#I^2 = 44.3% [0.0%; 83.4%]; H = 1.34 [1.00; 2.45]

#Test of heterogeneity:
#  Q d.f. p-value
#3.59    2  0.1662

###################
#First MR analysis#
###################

mr(dat_1_pre)

# id.exposure id.outcome outcome exposure                    method nsnp          b        se       pval
#1    moderate        BMI     BMI moderate                  MR Egger    3  0.1185481 0.4515826 0.83656431
#2    moderate        BMI     BMI moderate           Weighted median    3 -0.1771776 0.1175256 0.13166563
#3    moderate        BMI     BMI moderate Inverse variance weighted    3 -0.2363606 0.1225909 0.05384971
#4    moderate        BMI     BMI moderate               Simple mode    3 -0.1457169 0.1521639 0.43930615
#5    moderate        BMI     BMI moderate             Weighted mode    3 -0.1430262 0.1466260 0.43221610

#######################
#First Rucker Analysis#
#######################

mr_rucker(dat_1_pre)

#[[1]]$intercept
#Method    Estimate          SE      CI_low      CI_upp
#1  Egger fixed effects -0.01277599 0.009651194 -0.03169198 0.006140004
#2 Egger random effects -0.01277599 0.015525337 -0.04320509 0.017653113
#P
#1 0.1855785
#2 0.7947210

#[[1]]$Q
#Method        Q df         P
#1   Q_ivw 4.340111  2 0.1141713
#2 Q_egger 2.587736  1 0.1076942
#3  Q_diff 1.752375  1 0.1855785

#[[1]]$res
#[1] "A"

#[[1]]$selected
#Method nsnp   Estimate         SE    CI_low      CI_upp          P
#1 Rucker    3 -0.2363606 0.08321907 -0.399467 -0.07325422 0.00450822

#############
#CONCLUSIONS#
#############

dat_1_pre$labels <- NA

tiff("Moderate_BMI_Proxy_original_P07.png", units="in", width=10, height=10, res=300)
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

radial_output <- ivw_radial(radial_input, alpha = 0.05, weights = 3, tol = 0.0001) #NO OUTLIERS.

outliers <- radial_output$outliers$SNP #NO OUTLIERS.

#WE FINISH HERE.