#######
#CAUSE#
#######

#Inacc425alling CAUSE

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

###########
#FUNCTIONS#
###########

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

#########################
#Making space for memory#
#########################

memory.limit(size=8000000)

##############
#Loading data#
##############

WHRadjBMI <- fread("~/RAW_DATA/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
acc425 <- fread("~/RAW_DATA/Acc425_Model1_BOLTLMM_500K.txt.gz")

#########################
#Cleaning WHRadjBMI data#
#########################

summary(WHRadjBMI$Freq_Tested_Allele)

#MAF > 0.01
#INFO > 0.7
#Only keeping classic alleles (A, C, G and T)

WHRadjBMI <- WHRadjBMI[which(as.numeric(WHRadjBMI$Freq_Tested_Allele) > 0.01),]
WHRadjBMI <- WHRadjBMI[which(as.numeric(WHRadjBMI$Freq_Tested_Allele) < 0.99),]
WHRadjBMI <- WHRadjBMI[which(WHRadjBMI$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

WHRadjBMI <- WHRadjBMI[which(WHRadjBMI$Tested_Allele%in%yes_vect),]
WHRadjBMI <- WHRadjBMI[which(WHRadjBMI$Other_Allele%in%yes_vect),]

#chr6:28,477,797-33,448,354

WHRadjBMI_mhc <- WHRadjBMI[which(WHRadjBMI$CHR == 6 & WHRadjBMI$POS > 26000000 & WHRadjBMI$POS < 34000000),] #so we are gonna approximate it this way.

#This is the way that it is done by Warrington et al 2019.

WHRadjBMI_ <- WHRadjBMI[which(!(WHRadjBMI$SNP%in%WHRadjBMI_mhc$SNP)),]

WHRadjBMI <- WHRadjBMI_

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
WHRadjBMI$chr_pos <- paste(WHRadjBMI$CHR, WHRadjBMI$POS, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(WHRadjBMI, acc425, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("Tested_Allele", "ALLELE1"), 
                A2_cols = c("Other_Allele", "ALLELE0"))


#Checking that X is cool:

head(X)

acc425[which(acc425$chr_pos == "1:693731"),]
WHRadjBMI[which(WHRadjBMI$chr_pos == "1:693731"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

WHRadjBMI_match <- WHRadjBMI[which(WHRadjBMI$chr_pos%in%X$snp),]

WHRadjBMI_match <- WHRadjBMI_match[order(match(WHRadjBMI_match$chr_pos, X$snp)),]

which(WHRadjBMI_match$chr_pos != X$snp) #they all match. 

#Let's get the clean RSIDs:

WHRadjBMI_match$RSID <- as.character(unlist(sapply(WHRadjBMI_match$SNP, parse_rsid)))

X$snp <- WHRadjBMI_match$RSID

#Checked first one and last one in PhenoScanner to be sure.
#All good.

#Calculating Nuisance:

set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)

saveRDS(X, file = "~/X_WHRadjBMI_acc425.rds")
saveRDS(params, file = "~/params_WHRadjBMI_acc425.rds")

params <- readRDS("~/params_WHRadjBMI_acc425.rds")
X <- readRDS("~/X_WHRadjBMI_acc425.rds")

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

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

write.table(pruned_all, "~/WHRadjBMI_pruned_acc425.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

top_WHRadjBMI_pruned_df <- read.table("~/WHRadjBMI_pruned_acc425.txt", header = TRUE, stringsAsFactors = FALSE)
top_WHRadjBMI_pruned_snps <- top_WHRadjBMI_pruned_df$x

#Let's run it!!

res <- cause(X=X, variants = top_WHRadjBMI_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_WHRadjBMI_to_acc425.rds")

res <- readRDS("~/RES_WHRadjBMI_to_acc425_FULL.rds")

res$elpd

#model1  model2 delta_elpd se_delta_elpd          z
#1    null sharing 0.05689424     0.8383037 0.06786829
#2    null  causal 0.31263951     1.5205144 0.20561431
#3 sharing  causal 0.25574528     0.6829269 0.37448410

#There does not seem to be a causality...

summary(res)

#p-value testing that causal model is a better fit:  0.65 
#Posterior medians and  95 % credible intervals:
#  model     gamma                 eta                  q               
#[1,] "Sharing" NA                    "-0.04 (-0.3, 0.12)" "0.16 (0, 0.76)"
#[2,] "Causal"  "-0.01 (-0.06, 0.03)" "0 (-0.25, 0.21)"    "0.18 (0, 0.85)"

#################acc425 to WHRadjBMI################

#########################
#Making space for memory#
#########################

memory.limit(size=8000000)

##############
#Loading data#
##############

WHRadjBMI <- fread("~/RAW_DATA/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
acc425 <- fread("~/RAW_DATA/Acc425_Model1_BOLTLMM_500K.txt.gz")

#########################
#Cleaning WHRadjBMI data#
#########################

summary(WHRadjBMI$Freq_Tested_Allele)

#MAF > 0.01
#INFO > 0.7
#Only keeping classic alleles (A, C, G and T)

WHRadjBMI <- WHRadjBMI[which(as.numeric(WHRadjBMI$Freq_Tested_Allele) > 0.01),]
WHRadjBMI <- WHRadjBMI[which(as.numeric(WHRadjBMI$Freq_Tested_Allele) < 0.99),]
WHRadjBMI <- WHRadjBMI[which(WHRadjBMI$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

WHRadjBMI <- WHRadjBMI[which(WHRadjBMI$Tested_Allele%in%yes_vect),]
WHRadjBMI <- WHRadjBMI[which(WHRadjBMI$Other_Allele%in%yes_vect),]

#chr6:28,477,797-33,448,354

WHRadjBMI_mhc <- WHRadjBMI[which(WHRadjBMI$CHR == 6 & WHRadjBMI$POS > 26000000 & WHRadjBMI$POS < 34000000),] #so we are gonna approximate it this way.

#This is the way that it is done by Warrington et al 2019.

WHRadjBMI_ <- WHRadjBMI[which(!(WHRadjBMI$SNP%in%WHRadjBMI_mhc$SNP)),]

WHRadjBMI <- WHRadjBMI_

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
WHRadjBMI$chr_pos <- paste(WHRadjBMI$CHR, WHRadjBMI$POS, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(acc425, WHRadjBMI, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("ALLELE1", "Tested_Allele"), 
                A2_cols = c("ALLELE0", "Other_Allele"))


#Checking that X is cool:

head(X)

acc425[which(acc425$chr_pos == "1:693731"),]
WHRadjBMI[which(WHRadjBMI$chr_pos == "1:693731"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

WHRadjBMI_match <- WHRadjBMI[which(WHRadjBMI$chr_pos%in%X$snp),]

WHRadjBMI_match <- WHRadjBMI_match[order(match(WHRadjBMI_match$chr_pos, X$snp)),]

which(WHRadjBMI_match$chr_pos != X$snp) #they all match. 

head(WHRadjBMI_match$chr_pos)
head(X$snp)

#Let's get the clean RSIDs:

WHRadjBMI_match$RSID <- as.character(unlist(sapply(WHRadjBMI_match$SNP, parse_rsid)))

X$snp <- WHRadjBMI_match$RSID

#We got everything!!
#Let's save this X:

set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)

saveRDS(X, file = "~/X_acc425_WHRadjBMI.rds")
saveRDS(params, file = "~/params_acc425_WHRadjBMI.rds")

#check_X <- readRDS("~/LD_panels/LDShrink/CAUSE_Emergency_code/X_WHRadjBMI_acc425.rds")
#They are not gonna be the same due to the duplicates in the exposure, which changes per trait.
#So far so good.

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

#Since we used the RSIDs of WHRadjBMI, we are going to match the results with the WHRadjBMI dataframe.

WHRadjBMI_pruned <- WHRadjBMI_match[which(WHRadjBMI_match$RSID%in%pruned_all),]

acc425_pruned <- acc425[which(acc425$chr_pos%in%WHRadjBMI_pruned$chr_pos),]

write.table(acc425_pruned, "~/acc425_2_WHRadjBMI_pruned.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Now let's run the analysis!!!

top_acc425_pruned_df <- read.table("~/acc425_2_WHRadjBMI_pruned.txt", header = TRUE, stringsAsFactors = FALSE)

top_acc425_pruned_snps <- WHRadjBMI_match$RSID[which(WHRadjBMI_match$chr_pos%in%top_acc425_pruned_df$chr_pos)]

#Let's calcualte this!

res <- cause(X=X, variants = top_acc425_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_acc425_to_WHRadjBMI.rds")

#There seem to be problems. Let's check:

res$elpd

#model1  model2 delta_elpd se_delta_elpd          z
#1    null sharing -0.9142619     1.1899184 -0.7683400
#2    null  causal -1.6812004     2.0952777 -0.8023759
#3 sharing  causal -0.7669385     0.9075441 -0.8450702

summary(res)

#p-value testing that causal model is a better fit:  0.2 
#Posterior medians and  95 % credible intervals:
#  model     gamma                 eta                  q                  
#[1,] "Sharing" NA                    "-0.1 (-0.39, 0.07)" "0.36 (0.01, 0.86)"
#[2,] "Causal"  "-0.06 (-0.16, 0.04)" "0 (-0.37, 0.44)"    "0.19 (0, 0.86)"   

#plot(res, type="data")
