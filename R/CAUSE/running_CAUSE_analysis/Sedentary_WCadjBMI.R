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

##############
#Loading data#
##############

WCAdjBMI <- fread("~/CURATED_DATA/WCAdjBMI_combined_Curated.txt")
ST <- fread("~/RAW_DATA/Doherty-2018-NatureComms-sedentary.csv.gz")

####################
#Filtering WCAdjBMI#
####################

#We don't need to filter much because we have already done it before!

#We do not have info and I checked the frequencies. All good:

summary(WCAdjBMI$FreqAllele1HapMapCEU) #perfect.

#The main issue is gonna be... the matching.

##################
#Finally, alleles#
##################

yes_vect <- c("A", "G", "T", "C")

WCAdjBMI <- WCAdjBMI[which(WCAdjBMI$Allele1%in%yes_vect),] #all of them, as expected.
WCAdjBMI <- WCAdjBMI[which(WCAdjBMI$Allele2%in%yes_vect),] #all of them, as expected.

###############################################################################
#Let's do first a match of those SNPs that do not have chromosome and position#
###############################################################################

WCAdjBMI_missing <- WCAdjBMI[which(WCAdjBMI$chr_pos_37 == "-"),]

which(WCAdjBMI_missing$MarkerName%in%ST$SNP) #none of them. 

#Let's check the merged allele versions of those RSIDs:

merged_rs <- fread("~/RAW_DATA/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

###################
#Now we match them#
###################

index_merged_1 <- which(merged_rs$old_rs%in%WCAdjBMI_missing$MarkerName) #Only 2! 
index_merged_2 <- which(merged_rs$new_rs%in%WCAdjBMI_missing$MarkerName) #613!

swapping_snps_1 <- merged_rs[index_merged_1,]
swapping_snps_2 <- merged_rs[index_merged_2,]

which(swapping_snps_1$new_rs%in%ST$SNP) #0
which(swapping_snps_1$old_rs%in%ST$SNP) #0

which(swapping_snps_2$new_rs%in%ST$SNP) #0
which(swapping_snps_2$old_rs%in%ST$SNP) #0

#It doesn't matter which combination we do, this still gets out of luck.
#Hence, we are just gonna work with what we have: the chromosome and positions.

##################
#Cleaning ST data#
##################

#MAF > 0.01
#INFO > 0.7
#Also, removthing those that are deletions or insertions (classic alleles.)

summary(ST$A1FREQ) #maf is already done!
summary(ST$INFO) #INFO needs to be done:

ST<- ST[which(ST$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

ST<- ST[which(ST$ALLELE1%in%yes_vect),]
ST<- ST[which(ST$ALLELE0%in%yes_vect),]

###########################
#Removing MHC region in ST#
###########################

ST_mhc <- ST[which(ST$CHR == 6 & ST$BP >= 26000000 & ST$BP <= 34000000),]

summary(ST_mhc$CHR) #perfect
summary(ST_mhc$BP) #perfect

ST_ <- ST[-which(ST$SNP%in%ST_mhc$SNP),] #removes all the annoying ones perfectly.

ST <- ST_

####################################
#Matching data in the best possible#
####################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure firSTand do the matching later to save some of them.

ST$chr_pos <- paste("chr", ST$CHR, ":", ST$BP, sep = "")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(WCAdjBMI, ST, snp_name_cols = c("chr_pos_37", "chr_pos"), 
                beta_hat_cols = c("b", "BETA"), 
                se_cols = c("se", "SE"), 
                A1_cols = c("Allele1", "ALLELE1"), 
                A2_cols = c("Allele2", "ALLELE0"))


#Checking that X is cool:

head(X)

ST[which(ST$chr_pos == "chr7:92383888"),]
WCAdjBMI[which(WCAdjBMI$chr_pos_37 == "chr7:92383888"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

WCAdjBMI_match <- WCAdjBMI[which(WCAdjBMI$chr_pos_37%in%X$snp),]

WCAdjBMI_match <- WCAdjBMI_match[order(match(WCAdjBMI_match$chr_pos_37, X$snp)),]

which(WCAdjBMI_match$chr_pos_37 != X$snp) #they all match. 
length(which(WCAdjBMI_match$chr_pos_37 == X$snp)) #they all match. 

#Let's see if we have any weird IDs left:

check <- WCAdjBMI_match[order(WCAdjBMI_match$MarkerName),]

head(check, 300) #it seems like there are no mismatches!!
tail(check, 300) #it seems like there are no mismatches!!

index_clean <- which(str_detect(WCAdjBMI_match$MarkerName, "rs") != FALSE) #exactly!

#Let's follow the protocol with X_clean, so we can obtain the rsid
#and double check later.

X_clean <- X[index_clean,]
WCAdjBMI_clean <- WCAdjBMI_match[index_clean,]

which(X_clean$snp != WCAdjBMI_clean$chr_pos_37)

head(X_clean$snp)
head(WCAdjBMI_clean$chr_pos_37)

tail(X_clean$snp)
tail(WCAdjBMI_clean$chr_pos_37)

X_clean$snp <- WCAdjBMI_clean$MarkerName #now we can return the RSID.

#Let's change the name to keep up with the nomenclature of the other scripts:

X_end <- X_clean

head(X_end)
tail(X_end)

#Calculating Nuisance:

set.seed(100)
varlist<- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_WCAdjBMI_ST.rds")
saveRDS(params, file = "~/params_WCAdjBMI_ST.rds")

#params has been calculated, so let's go for LD pruning if possible:

library(tidyverse)

variants <- X_end %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

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

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

WCAdjBMI_pruned <- WCAdjBMI_match[which(WCAdjBMI_match$MarkerName%in%pruned_all),]

write.table(WCAdjBMI_pruned, "~/WCAdjBMI_pruned_ST.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

######################################################

top_WCAdjBMI_pruned_df <- read.table("~/WCAdjBMI_pruned_ST.txt", header = TRUE, stringsAsFactors = FALSE)
top_WCAdjBMI_pruned_snps <- top_WCAdjBMI_pruned_df$MarkerName

#Let's run it!!

res <- cause(X=X_end, variants = top_WCAdjBMI_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_WCAdjBMI_to_ST.rds")

res <- readRDS("~/RES_WCAdjBMI_to_ST.rds")

res$elpd

#  model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing  -7.588635     3.3439978 -2.269330
#2    null  causal  -8.608108     3.9697708 -2.168414
#3 sharing  causal  -1.019473     0.6536404 -1.559686

summary(res)

#p-value testing that causal model is a better fit:  0.059 
#Posterior medians and  95 % credible intervals:
#  model     gamma               eta                 q                  
#[1,] "Sharing" NA                  "0.1 (0.05, 0.2)"   "0.63 (0.23, 0.93)"
#[2,] "Causal"  "0.07 (0.01, 0.14)" "0.01 (-0.3, 0.48)" "0.17 (0, 0.85)"  

plot(res, type="data")

################################################################################
################################################################################
################################################################################

#################ST WCAdjBMI################

#########################
#Making space for memory#
#########################

memory.limit(size=80000000)

##############
#Loading data#
##############

WCAdjBMI <- fread("~/CURATED_DATA/WCAdjBMI_combined_Curated.txt")
ST <- fread("~/RAW_DATA/Doherty_Traits/Doherty-2018-NatureComms-sedentary.csv.gz")

####################
#Filtering WCAdjBMI#
####################

#We don't need to filter much because we have already done it before!

#We do not have info and I checked the frequencies. All good:

summary(WCAdjBMI$FreqAllele1HapMapCEU) #perfect.

#The main issue is gonna be... the matching.

##################
#Finally, alleles#
##################

yes_vect <- c("A", "G", "T", "C")

WCAdjBMI <- WCAdjBMI[which(WCAdjBMI$Allele1%in%yes_vect),] #all of them, as expected.
WCAdjBMI <- WCAdjBMI[which(WCAdjBMI$Allele2%in%yes_vect),] #all of them, as expected.

###############################################################################
#Let's do first a match of those SNPs that do not have chromosome and position#
###############################################################################

WCAdjBMI_missing <- WCAdjBMI[which(WCAdjBMI$chr_pos_37 == "-"),]

which(WCAdjBMI_missing$MarkerName%in%ST$SNP) #none of them. 

#Let's check the merged allele versions of those RSIDs:

merged_rs <- fread("~/RAW_DATA/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

###################
#Now we match them#
###################

index_merged_1 <- which(merged_rs$old_rs%in%WCAdjBMI_missing$MarkerName) #Only 2! 
index_merged_2 <- which(merged_rs$new_rs%in%WCAdjBMI_missing$MarkerName) #613!

swapping_snps_1 <- merged_rs[index_merged_1,]
swapping_snps_2 <- merged_rs[index_merged_2,]

which(swapping_snps_1$new_rs%in%ST$SNP) #0
which(swapping_snps_1$old_rs%in%ST$SNP) #0

which(swapping_snps_2$new_rs%in%ST$SNP) #0
which(swapping_snps_2$old_rs%in%ST$SNP) #0

#It doesn't matter which combination we do, this still gets out of luck.
#Hence, we are just gonna work with what we have: the chromosome and positions.

########################
#Cleaning ST data#
########################

#MAF > 0.01
#INFO > 0.7
#Also, removthing those that are deletions or insertions (classic alleles.)

summary(ST$A1FREQ) #maf is already done!
summary(ST$INFO) #INFO needs to be done:

ST<- ST[which(ST$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

ST<- ST[which(ST$ALLELE1%in%yes_vect),]
ST<- ST[which(ST$ALLELE0%in%yes_vect),]

###############################
#Removing MHC region in ST#
###############################

ST_mhc <- ST[which(ST$CHR == 6 & ST$BP >= 26000000 & ST$BP <= 34000000),]

summary(ST_mhc$CHR) #perfect
summary(ST_mhc$BP) #perfect

ST_ <- ST[-which(ST$SNP%in%ST_mhc$SNP),] #removes all the annoying ones perfectly.

ST <- ST_

####################################
#Matching data in the best possible#
####################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure firSTand do the matching later to save some of them.

ST$chr_pos <- paste("chr", ST$CHR, ":", ST$BP, sep = "")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(ST, WCAdjBMI, snp_name_cols = c("chr_pos", "chr_pos_37"), 
                beta_hat_cols = c("BETA", "b"), 
                se_cols = c("SE", "se"), 
                A1_cols = c("ALLELE1", "Allele1"), 
                A2_cols = c("ALLELE0", "Allele2"))


#Checking that X is cool:

head(X)

ST[which(ST$chr_pos == "chr1:752566"),]
WCAdjBMI[which(WCAdjBMI$chr_pos_37 == "chr1:752566"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

WCAdjBMI_match <- WCAdjBMI[which(WCAdjBMI$chr_pos_37%in%X$snp),]

WCAdjBMI_match <- WCAdjBMI_match[order(match(WCAdjBMI_match$chr_pos_37, X$snp)),]

which(WCAdjBMI_match$chr_pos_37 != X$snp) #they all match. 
length(which(WCAdjBMI_match$chr_pos_37 == X$snp)) #they all match. 

#Let's see if we have any weird IDs left:

check <- WCAdjBMI_match[order(WCAdjBMI_match$MarkerName),]

head(check, 300) #it seems like there are no mismatches!!
tail(check, 300) #it seems like there are no mismatches!!

index_clean <- which(str_detect(WCAdjBMI_match$MarkerName, "rs") != FALSE) #exactly!

#Let's follow the protocol with X_clean, so we can obtain the rsid
#and double check later.

X_clean <- X[index_clean,]
WCAdjBMI_clean <- WCAdjBMI_match[index_clean,]

which(X_clean$snp != WCAdjBMI_clean$chr_pos_37)

head(X_clean$snp)
head(WCAdjBMI_clean$chr_pos_37)

tail(X_clean$snp)
tail(WCAdjBMI_clean$chr_pos_37)

X_clean$snp <- WCAdjBMI_clean$MarkerName #now we can return the RSID.

#Let's change the name to keep up with the nomenclature of the other scripts:

X_end <- X_clean

head(X_end)
tail(X_end)

#Calculating Nuisance:

set.seed(100)
varlist<- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_ST_WCAdjBMI.rds")
saveRDS(params, file = "~/params_ST_WCAdjBMI.rds")

#Params has been calculated, so let's go for LD pruning if possible:

variants <- X_end %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

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

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

WCAdjBMI_pruned <- WCAdjBMI_match[which(WCAdjBMI_match$MarkerName%in%pruned_all),]

ST_pruned <- ST[which(ST$chr_pos%in%WCAdjBMI_pruned$chr_pos),]

write.table(ST_pruned, "~/ST_WCAdjBMI_pruned.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Let's retrieve the SNPs:

top_ST_pruned_df <- fread("~/ST_WCAdjBMI_pruned.txt")

top_ST_pruned_snps <- top_ST_pruned_df$SNP

#Let's calcualte this!

res <- cause(X=X_end, variants = top_ST_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_ST_to_WCAdjBMI.rds")

res <- readRDS("~/RES_ST_to_WCAdjBMI.rds")

#There seem to be problems. Let's check:

res$elpd

#model1  model2 delta_elpd se_delta_elpd        z
#1    null sharing  0.5549051    0.07338547 7.561512
#2    null  causal  1.3892542    0.16248979 8.549794
#3 sharing  causal  0.8343491    0.09060456 9.208688

summary(res)

#p-value testing that causal model is a better fit:  1 
#Posterior medians and  95 % credible intervals:
#  model     gamma             eta                  q                  
#[1,] "Sharing" NA                "0.01 (-0.44, 0.57)" "0.12 (0, 0.71)"   
#[2,] "Causal"  "0 (-0.13, 0.14)" "0.01 (-0.41, 0.52)" "0.21 (0.01, 0.85)"

plot(res, type="data")