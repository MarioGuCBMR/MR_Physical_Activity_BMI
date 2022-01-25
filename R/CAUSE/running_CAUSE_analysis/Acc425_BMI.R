#######
#CAUSE#
#######

#Installing CAUSE

#inacc425all.packages("devtools")
#inacc425all.packages("tidyverse")
#devtools::inacc425all_version("mixsqp", version = "0.1-97", repos = "http://cran.us.r-project.org")
#devtools::inacc425all_version("ashr", version = "2.2-32", repos = "http://cran.us.r-project.org")
#devtools::inacc425all_github("jean997/cause")

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

memory.limit(size=8000000)

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
acc425 <- fread("~/RAW_DATA/Acc425_Model1_BOLTLMM_500K.txt.gz")

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

######################
#Cleaning acc425 data#
######################

#MAF > 0.01
#INFO > 0.7
#Also, removthing those that are deletions or insertions (classic alleles.)

summary(acc425$A1FREQ) #For acc425 we have to change it!!

acc425 <- acc425[which(acc425$A1FREQ > 0.01),]
acc425 <- acc425[which(acc425$A1FREQ < 0.99),]

summary(acc425$INFO) #INFO needs to be done:

acc425 <- acc425[which(acc425$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

acc425 <- acc425[which(acc425$ALLELE1%in%yes_vect),]
acc425 <- acc425[which(acc425$ALLELE0%in%yes_vect),]

########################################
#Matching data in the best way possible#
########################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure firacc425 and do the matching later to save some of them.

acc425$chr_pos <- paste(acc425$CHR, acc425$BP, sep = ":")
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

X <- gwas_merge(BMI, acc425, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("A1", "ALLELE1"), 
                A2_cols = c("A2", "ALLELE0"))


#Checking that X is cool:

head(X)

acc425[which(acc425$chr_pos == "1:693731"),]
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

saveRDS(X, file = "~/X_BMI_acc425.rds")
saveRDS(params, file = "~/params_BMI_acc425.rds")

params <- readRDS("~/params_BMI_acc425.rds")
X <- readRDS("~/X_BMI_acc425.rds")

#params has been calculated, so let's go for LD pruning if possible:

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

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", deacc425file = "~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", deacc425file="~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

BMI_pruned <- BMI_match[which(BMI_match$SNP%in%pruned_all),]

write.table(BMI_pruned, "~/BMI_pruned_acc425.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Let's retrieve the SNPs:

top_bmi_pruned_df <- read.table("~/BMI_pruned_acc425.txt", header = TRUE, stringsAsFactors = FALSE)
top_bmi_pruned_snps <- top_bmi_pruned_df$SNP

#Let's run it!!

res <- cause(X=X, variants = top_bmi_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_BMI_to_acc425.rds")

res <- readRDS("~/RES_BMI_to_acc425.rds")

res$elpd

#model1  model2   delta_elpd se_delta_elpd          z
#1    null sharing -251.8775693    21.5435795 -11.691538
#2    null  causal -252.2229798    21.6523330 -11.648767
#3 sharing  causal   -0.3454106     0.9018152  -0.383017

#There does not seem to be a causality...

summary(res)

#p-value testing that causal model is a better fit:  0.35 
#Posterior medians and  95 % credible intervals:
#  model     gamma                 eta                    q                  
#[1,] "Sharing" NA                    "-0.16 (-0.19, -0.14)" "0.9 (0.77, 0.98)" 
#[2,] "Causal"  "-0.15 (-0.2, -0.08)" "0.01 (-0.18, 0.21)"   "0.21 (0.01, 0.83)"

#################acc425 to BMI################

#########################
#Making space for memory#
#########################

memory.limit(size=8000000)

##############
#Loading data#
##############

BMI <- fread("~/RAW_DATA/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
acc425 <- fread("~/RAW_DATA/Acc425_Model1_BOLTLMM_500K.txt.gz")

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

######################
#Cleaning acc425 data#
######################

#MAF > 0.01
#INFO > 0.7
#Also, removthing those that are deletions or insertions (classic alleles.)

summary(acc425$A1FREQ) #maf is already done!
summary(acc425$INFO) #INFO needs to be done:

acc425 <- acc425[which(acc425$A1FREQ > 0.01),]
acc425 <- acc425[which(acc425$A1FREQ < 0.99),]
acc425 <- acc425[which(acc425$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

acc425 <- acc425[which(acc425$ALLELE1%in%yes_vect),]
acc425 <- acc425[which(acc425$ALLELE0%in%yes_vect),]

########################################
#Matching data in the best way possible#
########################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure firacc425 and do the matching later to save some of them.

acc425$chr_pos <- paste(acc425$CHR, acc425$BP, sep = ":")
BMI$chr_pos <- paste(BMI$CHR, BMI$POS, sep = ":")

#Getting the alleles in the increasing order for acc425, just in case.

new_a1 <- ifelse(acc425$BETA < 0, toupper(as.character(acc425$ALLELE0)), toupper(as.character(acc425$ALLELE1)))
new_a2 <- ifelse(acc425$BETA < 0, toupper(as.character(acc425$ALLELE1)), toupper(as.character(acc425$ALLELE0)))
beta <- ifelse(acc425$BETA < 0, acc425$BETA*(-1), acc425$BETA)

acc425$A1 <- new_a1
acc425$A0 <- new_a2
acc425$BETA_ <- beta

#acc425$ALLELE1 <- new_a1
#acc425$ALLELE0 <- new_a2
#acc425$BETA <- beta

#Careful there is some shit overhere

X <- gwas_merge(acc425, BMI, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA_", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("A1", "Tested_Allele"), 
                A2_cols = c("A0", "Other_Allele"))

#############
#Small check#
#############

head(X)
BMI[which(BMI$chr_pos == "1:693731")]
acc425[which(acc425$chr_pos == "1:693731")]

#Now we use BMI to get the RSIDs:

BMI_match <- BMI[which(BMI$chr_pos%in%X$snp),]

BMI_match <- BMI_match[order(match(BMI_match$chr_pos, X$snp)),]

#Small check:

which(BMI_match$chr_pos != X$snp) #they all match:

BMI_match$SNP <- as.character(unlist(sapply(BMI_match$SNP, parse_rsid)))

X$snp <- BMI_match$SNP

#Double-checked. All the conversion are done properly and the switching of the strands and all is done correctly.

#We got everything!!
#Let's save this X:

set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)

saveRDS(X, file = "~/X_acc425_BMI.rds")
saveRDS(params, file = "~/params_acc425_BMI.rds")

params_ <- readRDS("~/params_acc425_BMI.rds")
X_ <- readRDS("~/X_acc425_BMI.rds")

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

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", deacc425file = "~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", deacc425file="~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

#Careful, we lose one. 

BMI_pruned <- BMI_match[which(BMI_match$SNP%in%pruned_all),]

acc425_pruned <- acc425[which(acc425$chr_pos%in%BMI_pruned$chr_pos),]

write.table(acc425_pruned, "~/acc425_pruned.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Let's get the SNPs:

top_acc425_pruned_df <- read.table("~/acc425_pruned.txt", header = TRUE, stringsAsFactors = FALSE)

#Just in case we might miss the weird SNP that we also detect in Moderate, we are going to use BMI_match and its RSIDs:

top_acc425_pruned_snps <- BMI_match$SNP[which(BMI_match$chr_pos%in%top_acc425_pruned_df$chr_pos)]

#Let's calcualte this!

res <- cause(X=X, variants = top_acc425_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_acc425_to_BMI.rds")

#There seem to be problems. Let's check:

res$elpd

#model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing -17.472808     4.2397015 -4.121235
#2    null  causal -18.981530     4.6821749 -4.053999
#3 sharing  causal  -1.508722     0.4484551 -3.364265

summary(res)

#p-value testing that causal model is a better fit:  0.00038 
#Posterior medians and  95 % credible intervals:
#  model     gamma                  eta                    q                  
#[1,] "Sharing" NA                     "-0.18 (-0.25, -0.12)" "0.85 (0.57, 0.98)"
#[2,] "Causal"  "-0.16 (-0.24, -0.08)" "0 (-0.31, 0.31)"      "0.19 (0, 0.86)"   

#plot(res, type="data")
