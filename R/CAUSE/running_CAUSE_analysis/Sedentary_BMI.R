#######
#CAUSE#
#######

#Installing CAUSE

#install.packages("devtools")
#install.packages("tidyverse")
#devtools::install_version("mixsqp", version = "0.1-97", repos = "http://cran.us.r-project.org")
#devtools::install_version("ashr", version = "2.2-32", repos = "http://cran.us.r-project.org")
#devtools::install_github("jean997/cause")

#################
#Loading libries#
#################

library(data.table)
library(cause)
library(tidyverse)
library(TwoSampleMR)

#########################
#Making space for memory#
#########################

memory.limit(size=80000000)

###################
#Loading functions#
###################

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

##############
#Loading data#
##############

BMI <- fread("~/RAW_DATA/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
ST <- fread("~/RAW_DATA/Doherty_Traits/Doherty-2018-NatureComms-sedentary.csv.gz")

###################
#Cleaning BMI data#
###################

summary(BMI$Freq_Tested_Allele)

#MAF > 0.01
#INFO > 0.7
#Only keeping classic alleles (A, C, G and T)

BMI <- BMI[which(as.numeric(BMI$Freq_Tested_Allele) > 0.01),]
BMI <- BMI[which(as.numeric(BMI$Freq_Tested_Allele) < 0.99),]
BMI <- BMI[which(BMI$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

BMI <- BMI[which(BMI$Tested_Allele%in%yes_vect),]
BMI <- BMI[which(BMI$Other_Allele%in%yes_vect),]

#And now we remove the data in the MHC region:

#The expanded MHC is 26M to 24M.

BMI_mhc <- BMI[which(BMI$CHR == 6 & BMI$POS > 26000000 & BMI$POS < 34000000),] #so we are gonna approximate it this way.

#Small  check:

summary(BMI_mhc$CHR) #perfect
summary(BMI_mhc$POS) #perfect

#This is the way that it is done by Warrington et al 2019.

BMI_ <- BMI[which(!(BMI$SNP%in%BMI_mhc$SNP)),]

BMI <- BMI_

#Perfect!

#########################
#Cleaning Sedentary data#
#########################

#MAF > 0.01
#INFO > 0.7
#Also, removthing those that are deletions or insertions (classic alleles.)

summary(ST$A1FREQ) #maf is already done!
summary(ST$INFO) #INFO needs to be done:

ST <- ST[which(ST$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

ST <- ST[which(ST$ALLELE1%in%yes_vect),]
ST <- ST[which(ST$ALLELE0%in%yes_vect),]

########################################
#Matching data in the best way possible#
########################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure first and do the matching later to save some of them.

ST$chr_pos <- paste(ST$CHR, ST$BP, sep = ":")
BMI$chr_pos <- paste(BMI$CHR, BMI$POS, sep = ":")

#######
#CAUSE#
#######

library(cause)

new_a1 <- ifelse(BMI$BETA < 0, toupper(as.character(BMI$Other_Allele)), toupper(as.character(BMI$Tested_Allele)))
new_a2 <- ifelse(BMI$BETA < 0, toupper(as.character(BMI$Tested_Allele)), toupper(as.character(BMI$Other_Allele)))
beta <- ifelse(BMI$BETA < 0, BMI$BETA*(-1), BMI$BETA)

BMI$A1 <- new_a1
BMI$A2 <- new_a2
BMI$BETA <- beta

X <- gwas_merge(BMI, ST, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("A1", "ALLELE1"), 
                A2_cols = c("A2", "ALLELE0"))


#Checking that X is cool:

head(X)

ST[which(ST$chr_pos == "1:693731"),]
BMI[which(BMI$chr_pos == "1:693731"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

BMI_match <- BMI[which(BMI$chr_pos%in%X$snp),]

BMI_match <- BMI_match[order(match(BMI_match$chr_pos, X$snp)),]

which(BMI_match$chr_pos != X$snp) #they all match. 

#Checked first one and last one in PhenoScanner to be sure.
#All good.

#Let's clean the BMI RSIDs:

BMI_match$SNP <- as.character(unlist(sapply(BMI_match$SNP, parse_rsid)))

X$snp <- BMI_match$SNP

#Calculating Nuisance:

set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)

saveRDS(X, file = "~/X_BMI_ST.rds")
saveRDS(params, file = "~/params_BMI_ST_ORIGINAL.rds")

#STrams has been calculated, so let's go for LD pruning if possible:

library(tidyverse)

variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

#1

ld <- readRDS("~/LD_panels/LDShrink/chr1_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr1_AF0.05_snpdata.RDS")

pruned_1 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#2

ld <- readRDS("~/LD_panels/LDShrink/chr2_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr2_AF0.05_snpdata.RDS")

pruned_2 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#3

ld <- readRDS("~/LD_panels/LDShrink/chr3_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr3_AF0.05_snpdata.RDS")

pruned_3 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#4

ld <- readRDS("~/LD_panels/LDShrink/chr4_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr4_AF0.05_snpdata.RDS")

pruned_4 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#5

ld <- readRDS("~/LD_panels/LDShrink/chr5_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr5_AF0.05_snpdata.RDS")

pruned_5 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#6

ld <- readRDS("~/LD_panels/LDShrink/chr6_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr6_AF0.05_snpdata.RDS")

pruned_6 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#7

ld <- readRDS("~/LD_panels/LDShrink/chr7_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr7_AF0.05_snpdata.RDS")

pruned_7 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#8

ld <- readRDS("~/LD_panels/LDShrink/chr8_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr8_AF0.05_snpdata.RDS")

pruned_8 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#9

ld <- readRDS("~/LD_panels/LDShrink/chr9_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr9_AF0.05_snpdata.RDS")

pruned_9 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#10

ld <- readRDS("~/LD_panels/LDShrink/chr10_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr10_AF0.05_snpdata.RDS")

pruned_10 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#11

ld <- readRDS("~/LD_panels/LDShrink/chr11_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr11_AF0.05_snpdata.RDS")

pruned_11 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#12

ld <- readRDS("~/LD_panels/LDShrink/chr12_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr12_AF0.05_snpdata.RDS")

pruned_12 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#13

ld <- readRDS("~/LD_panels/LDShrink/chr13_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr13_AF0.05_snpdata.RDS")

pruned_13 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#14

ld <- readRDS("~/LD_panels/LDShrink/chr14_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr14_AF0.05_snpdata.RDS")

pruned_14 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#15

ld <- readRDS("~/LD_panels/LDShrink/chr15_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr15_AF0.05_snpdata.RDS")

pruned_15 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#16

ld <- readRDS("~/LD_panels/LDShrink/chr16_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr16_AF0.05_snpdata.RDS")

pruned_16 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#17

ld <- readRDS("~/LD_panels/LDShrink/chr17_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr17_AF0.05_snpdata.RDS")

pruned_17 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#18

ld <- readRDS("~/LD_panels/LDShrink/chr18_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr18_AF0.05_snpdata.RDS")

pruned_18 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#19

ld <- readRDS("~/LD_panels/LDShrink/chr19_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr19_AF0.05_snpdata.RDS")

pruned_19 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#20

ld <- readRDS("~/LD_panels/LDShrink/chr20_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr20_AF0.05_snpdata.RDS")

pruned_20 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#21

ld <- readRDS("~/LD_panels/LDShrink/chr21_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr21_AF0.05_snpdata.RDS")

pruned_21 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))
#22

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "~/RAW_DATA/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="~/RAW_DATA/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

BMI_pruned <- BMI_match[which(BMI_match$SNP%in%pruned_all),]

write.table(BMI_pruned, "~/BMI_pruned_ST.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

###############################################################################
#Let's load them again so that we can work with less memory after restarting R#
###############################################################################

top_bmi_pruned_df <- read.table("~/BMI_pruned_ST.txt", header = TRUE, stringsAsFactors = FALSE)
top_bmi_pruned_snps <- top_bmi_pruned_df$SNP

#Let's run it!!

res <- cause(X=X, variants = top_bmi_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_BMI_to_ST.rds")

res <- readRDS("~/RES_BMI_to_ST.rds")

res$elpd

#model1  model2  delta_elpd se_delta_elpd         z
#1    null sharing -137.838562    15.7463639 -8.753676
#2    null  causal -138.894977    15.8803011 -8.746369
#3 sharing  causal   -1.056415     0.3276895 -3.223828

summary(res)

#p-value testing that causal model is a better fit:  0.00063 
#Posterior medians and  95 % credible intervals:
#  model     gamma               eta                 q                 
#[1,] "Sharing" NA                  "0.14 (0.12, 0.17)" "0.9 (0.73, 0.98)"
#[2,] "Causal"  "0.13 (0.08, 0.17)" "0 (-0.23, 0.26)"   "0.18 (0, 0.85)"  

plot(res, type="data")

################################################################################
################################################################################
################################################################################

#################ST to BMI################

#########################
#Making space for memory#
#########################

memory.limit(size=800000)

##############
#Loading data#
##############

BMI <- fread("~/RAW_DATA/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
ST <- fread("~/RAW_DATA/Doherty_Traits/Doherty-2018-NatureComms-sedentary.csv.gz")

###################
#Cleaning BMI data#
###################

summary(BMI$Freq_Tested_Allele)

#MAF > 0.01
#INFO > 0.7
#Only keeping classic alleles (A, C, G and T)

BMI <- BMI[which(as.numeric(BMI$Freq_Tested_Allele) > 0.01),]
BMI <- BMI[which(as.numeric(BMI$Freq_Tested_Allele) < 0.99),]
BMI <- BMI[which(BMI$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

BMI <- BMI[which(BMI$Tested_Allele%in%yes_vect),]
BMI <- BMI[which(BMI$Other_Allele%in%yes_vect),]

#And now we remove the data in the MHC region:

#The expanded MHC is 26M to 24M.

BMI_mhc <- BMI[which(BMI$CHR == 6 & BMI$POS > 26000000 & BMI$POS < 34000000),] #so we are gonna approximate it this way.

#Small  check:

summary(BMI_mhc$CHR) #perfect
summary(BMI_mhc$POS) #perfect

#This is the way that it is done by Warrington et al 2019.

BMI_ <- BMI[which(!(BMI$SNP%in%BMI_mhc$SNP)),]

BMI <- BMI_

#########################
#Cleaning Sedentary data#
#########################

#MAF > 0.01
#INFO > 0.7
#Also, removthing those that are deletions or insertions (classic alleles.)

summary(ST$A1FREQ) #maf is already done!
summary(ST$INFO) #INFO needs to be done:

ST <- ST[which(ST$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

ST <- ST[which(ST$ALLELE1%in%yes_vect),]
ST <- ST[which(ST$ALLELE0%in%yes_vect),]

########################################
#Matching data in the best way possible#
########################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure first and do the matching later to save some of them.

ST$chr_pos <- paste(ST$CHR, ST$BP, sep = ":")
BMI$chr_pos <- paste(BMI$CHR, BMI$POS, sep = ":")

#Getting the alleles in the increasing order for ST, just in case.

new_a1 <- ifelse(ST$BETA < 0, toupper(as.character(ST$ALLELE0)), toupper(as.character(ST$ALLELE1)))
new_a2 <- ifelse(ST$BETA < 0, toupper(as.character(ST$ALLELE1)), toupper(as.character(ST$ALLELE0)))
beta <- ifelse(ST$BETA < 0, ST$BETA*(-1), ST$BETA)

ST$ALLELE1 <- new_a1
ST$ALLELE0 <- new_a2
ST$BETA <- beta

#Careful there is some shit overhere

X <- gwas_merge(ST, BMI, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("ALLELE1", "Tested_Allele"), 
                A2_cols = c("ALLELE0", "Other_Allele"))

#Now we use BMI to get the RSIDs:

BMI_match <- BMI[which(BMI$chr_pos%in%X$snp),]

BMI_match <- BMI_match[order(match(BMI_match$chr_pos, X$snp)),]

#Small check:

which(BMI_match$chr_pos != X$snp) #they all match:

BMI_match$SNP <- as.character(unlist(sapply(BMI_match$SNP, parse_rsid)))

X$snp <- BMI_match$SNP

#We got everything!!
#Let's save this X:

set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)

saveRDS(X, file = " ~/X_ST_BMI.rds")
saveRDS(params, file = "~/params_ST_BMI.rds")

#Params has been calculated, so let's go for LD pruning if possible:

variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

#1

ld <- readRDS("~/LD_panels/LDShrink/chr1_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr1_AF0.05_snpdata.RDS")

pruned_1 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#2

ld <- readRDS("~/LD_panels/LDShrink/chr2_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr2_AF0.05_snpdata.RDS")

pruned_2 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#3

ld <- readRDS("~/LD_panels/LDShrink/chr3_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr3_AF0.05_snpdata.RDS")

pruned_3 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#4

ld <- readRDS("~/LD_panels/LDShrink/chr4_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr4_AF0.05_snpdata.RDS")

pruned_4 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#5

ld <- readRDS("~/LD_panels/LDShrink/chr5_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr5_AF0.05_snpdata.RDS")

pruned_5 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#6

ld <- readRDS("~/LD_panels/LDShrink/chr6_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr6_AF0.05_snpdata.RDS")

pruned_6 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#7

ld <- readRDS("~/LD_panels/LDShrink/chr7_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr7_AF0.05_snpdata.RDS")

pruned_7 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#8

ld <- readRDS("~/LD_panels/LDShrink/chr8_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr8_AF0.05_snpdata.RDS")

pruned_8 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#9

ld <- readRDS("~/LD_panels/LDShrink/chr9_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr9_AF0.05_snpdata.RDS")

pruned_9 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#10

ld <- readRDS("~/LD_panels/LDShrink/chr10_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr10_AF0.05_snpdata.RDS")

pruned_10 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#11

ld <- readRDS("~/LD_panels/LDShrink/chr11_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr11_AF0.05_snpdata.RDS")

pruned_11 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#12

ld <- readRDS("~/LD_panels/LDShrink/chr12_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr12_AF0.05_snpdata.RDS")

pruned_12 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#13

ld <- readRDS("~/LD_panels/LDShrink/chr13_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr13_AF0.05_snpdata.RDS")

pruned_13 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#14

ld <- readRDS("~/LD_panels/LDShrink/chr14_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr14_AF0.05_snpdata.RDS")

pruned_14 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#15

ld <- readRDS("~/LD_panels/LDShrink/chr15_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr15_AF0.05_snpdata.RDS")

pruned_15 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#16

ld <- readRDS("~/LD_panels/LDShrink/chr16_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr16_AF0.05_snpdata.RDS")

pruned_16 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#17

ld <- readRDS("~/LD_panels/LDShrink/chr17_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr17_AF0.05_snpdata.RDS")

pruned_17 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#18

ld <- readRDS("~/LD_panels/LDShrink/chr18_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr18_AF0.05_snpdata.RDS")

pruned_18 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#19

ld <- readRDS("~/LD_panels/LDShrink/chr19_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr19_AF0.05_snpdata.RDS")

pruned_19 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#20

ld <- readRDS("~/LD_panels/LDShrink/chr20_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr20_AF0.05_snpdata.RDS")

pruned_20 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#21

ld <- readRDS("~/LD_panels/LDShrink/chr21_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr21_AF0.05_snpdata.RDS")

pruned_21 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))
#22

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "~/RAW_DATA/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="~/RAW_DATA/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

BMI_pruned <- BMI_match[which(BMI_match$SNP%in%pruned_all),]

ST_pruned <- ST[which(ST$chr_pos%in%BMI_pruned$chr_pos),]

write.table(ST_pruned, "~/ST_BMI_pruned.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

top_ST_pruned_df <- fread("~/ST_BMI_pruned.txt")

#Let's get the SNPs:

top_ST_pruned_snps <- top_ST_pruned_df$SNP

#Let's calcualte this!

res <- cause(X=X, variants = top_ST_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_ST_to_BMI.rds")

res <- readRDS("~/RES_ST_to_BMI.rds")

res$elpd

#model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing  -6.507066     2.5747790 -2.527233
#2    null  causal  -7.669581     3.1375277 -2.444466
#3 sharing  causal  -1.162515     0.5659824 -2.053977

summary(res)

#p-value testing that causal model is a better fit:  0.02 
#Posterior medians and  95 % credible intervals:
#  model     gamma              eta                 q                  
#[1,] "Sharing" NA                 "0.13 (0.06, 0.23)" "0.69 (0.25, 0.95)"
#[2,] "Causal"  "0.11 (0.02, 0.2)" "0 (-0.32, 0.32)"   "0.19 (0, 0.86)"   

plot(res, type="data")