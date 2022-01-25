#######
#CAUSE#
#######

#Installing CAUSE

#install.packages("devtools")
#install.packages("tidyverse")
#devtools::install_version("mixsqp", version = "0.1-97", repos = "http://cran.us.r-project.org")
#devtools::install_version("ashr", version = "2.2-32", repos = "http://cran.us.r-project.org")
#devtools::install_github("jean997/cause")

###################
#Loading libraries#
###################

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

WholeBodyFatFreeMass <- fread("~/CURATED_DATA/whole_fat_free_mass_UKBB_Combined_Curated_FULL.txt")
acc425 <- fread("~/RAW_DATA/Acc425_Model1_BOLTLMM_500K.txt.gz")

########################################
#Cleaning Whole Body Fat Free Mass data#
########################################

#(Since we had to transform the data from vcf to a readble txt with all the necessary column)
#all of the changes that we are about to make have already been done!)
#Nonetheless, I will include them either way.

summary(WholeBodyFatFreeMass$eaf.exposure) 

#We do not have INFO so we skip that..., but we do have alleles! CAUSE will only detect
#standard alleles

yes_vect <- c("A", "G", "C", "T")

WholeBodyFatFreeMass <- WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$effect_allele.exposure%in%yes_vect),]
WholeBodyFatFreeMass <- WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$other_allele.exposure%in%yes_vect),]

#And now we remove the data in the MHC region:

#The expanded MHC is 26M to 24M.

WholeBodyFatFreeMass_mhc <- WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$chr.exposure == 6 & WholeBodyFatFreeMass$pos.exposure > 26000000 & WholeBodyFatFreeMass$pos.exposure < 34000000),] #so we are gonna approximate it this way.

#All of them were removed in the curation!!

#Perfect!

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
acc425<- acc425[which(acc425$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

acc425<- acc425[which(acc425$ALLELE1%in%yes_vect),]
acc425<- acc425[which(acc425$ALLELE0%in%yes_vect),]

####################################
#Matching data in the best possible#
####################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure firacc425and do the matching later to save some of them.

acc425$chr_pos <- paste(acc425$CHR, acc425$BP, sep = ":")
WholeBodyFatFreeMass$chr_pos <- paste(WholeBodyFatFreeMass$chr.exposure, WholeBodyFatFreeMass$pos.exposure, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(WholeBodyFatFreeMass, acc425, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("beta.exposure", "BETA"), 
                se_cols = c("se.exposure", "SE"), 
                A1_cols = c("effect_allele.exposure", "ALLELE1"), 
                A2_cols = c("other_allele.exposure", "ALLELE0"))


#Checking that X is cool:

head(X)

acc425[which(acc425$chr_pos == "7:92383888"),]
WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$chr_pos == "7:92383888"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

WholeBodyFatFreeMass_match <- WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$chr_pos%in%X$snp),]

WholeBodyFatFreeMass_match <- WholeBodyFatFreeMass_match[order(match(WholeBodyFatFreeMass_match$chr_pos, X$snp)),]

which(WholeBodyFatFreeMass_match$chr_pos != X$snp) #they all match. 

#We have 530 variants that do not possess rsIDs.

check <- WholeBodyFatFreeMass_match[order(WholeBodyFatFreeMass_match$SNP),]

head(check, 300) #it seems like there are no mismatches!!

index_clean <- which(str_detect(WholeBodyFatFreeMass_match$SNP, "rs") != FALSE) #exactly!

#Let's follow the protocol with X_clean, so we can obtain the rsid
#and double check later.

X_clean <- X[index_clean,]
WholeBodyFatFreeMass_clean <- WholeBodyFatFreeMass_match[index_clean,]

which(X_clean$snp != WholeBodyFatFreeMass_clean$chr_pos)

head(X_clean$snp)
head(WholeBodyFatFreeMass_clean$chr_pos)

tail(X_clean$snp)
tail(WholeBodyFatFreeMass_clean$chr_pos)

X_clean$snp <- WholeBodyFatFreeMass_clean$SNP #now we can return the RSID.

#Let's change the name to keep up with the nomenclature of the other scripts:

X_end <- X_clean

head(X_end)
tail(X_end)

#Calculating Nuisance:

set.seed(100)
varlist<- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_WholeBodyFatFreeMass_acc425.rds")
saveRDS(params, file = "~/params_WholeBodyFatFreeMass_acc425.rds")

params <- readRDS("~/params_WholeBodyFatFreeMass_acc425.rds")
X_end <- readRDS("~/X_WholeBodyFatFreeMass_acc425.rds")

#acc425params has been calculated, so let's go for LD pruning if possible:

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

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

WholeBodyFatFreeMass_pruned <- WholeBodyFatFreeMass_match[which(WholeBodyFatFreeMass_match$SNP%in%pruned_all),]

write.table(WholeBodyFatFreeMass_pruned, "~/WholeBodyFatFreeMass_pruned_acc425.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

######################################################

top_WholeBodyFatFreeMass_pruned_df <- read.table("~/WholeBodyFatFreeMass_pruned_acc425.txt", header = TRUE, stringsAsFactors = FALSE)
top_WholeBodyFatFreeMass_pruned_snps <- top_WholeBodyFatFreeMass_pruned_df$SNP

#Let's run it!!

res <- cause(X=X_end, variants = top_WholeBodyFatFreeMass_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_WholeBodyFatFreeMass_to_acc425.rds")

res <- readRDS("~/RES_WholeBodyFatFreeMass_to_acc425.rds")

res$elpd

#model1  model2  delta_elpd se_delta_elpd         z
#1    null sharing -80.9534171    12.0871167 -6.697496
#2    null  causal -81.8379239    12.3320502 -6.636198
#3 sharing  causal  -0.8845068     0.3304377 -2.676773

summary(res)

#p-value testing that causal model is a better fit:  0.0037 
#Posterior medians and  95 % credible intervals:
#  model     gamma                  eta                    q                  
#[1,] "Sharing" NA                     "-0.13 (-0.18, -0.11)" "0.85 (0.62, 0.97)"
#[2,] "Causal"  "-0.11 (-0.17, -0.06)" "0 (-0.25, 0.27)"      "0.19 (0, 0.84)"   

plot(res, type="data")

################################################################################
################################################################################
################################################################################

#################acc425 WholeBodyFatFreeMass################

#########################
#Making space for memory#
#########################

memory.limit(size=800000)

##############
#Loading data#
##############

WholeBodyFatFreeMass <- fread("~/CURATED_DATA/whole_fat_free_mass_UKBB_Combined_Curated_FULL.txt")
acc425 <- fread("~/RAW_DATA/Acc425_Model1_BOLTLMM_500K.txt.gz")

####################################
#Cleaning WholeBodyFatFreeMass data#
####################################

#(Since we had to transform the data from vcf to a readble txt with all the necessary column)
#all of the changes that we are about to make have already been done!)
#Nonetheless, I will include them either way.

summary(WholeBodyFatFreeMass$eaf.exposure) 

#We do not have INFO so we skip that..., but we do have alleles! CAUSE will only detect
#standard alleles

yes_vect <- c("A", "G", "C", "T")

WholeBodyFatFreeMass <- WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$effect_allele.exposure%in%yes_vect),]
WholeBodyFatFreeMass <- WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$other_allele.exposure%in%yes_vect),]

#And now we remove the data in the MHC region:

#The expanded MHC is 26M to 24M.

WholeBodyFatFreeMass_mhc <- WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$chr.exposure == 6 & WholeBodyFatFreeMass$pos.exposure > 26000000 & WholeBodyFatFreeMass$pos.exposure < 34000000),] #so we are gonna approximate it this way.

#All of them were removed in the curation!!

#Perfect!

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
acc425<- acc425[which(acc425$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

acc425<- acc425[which(acc425$ALLELE1%in%yes_vect),]
acc425<- acc425[which(acc425$ALLELE0%in%yes_vect),]


########################################
#Matching data in the best way possible#
########################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure first and do the matching later to save some of them.

acc425$chr_pos <- paste(acc425$CHR, acc425$BP, sep = ":")
WholeBodyFatFreeMass$chr_pos <- paste(WholeBodyFatFreeMass$chr.exposure, WholeBodyFatFreeMass$pos.exposure, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(acc425, WholeBodyFatFreeMass, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "beta.exposure"), 
                se_cols = c("SE", "se.exposure"), 
                A1_cols = c("ALLELE1", "effect_allele.exposure"), 
                A2_cols = c("ALLELE0", "other_allele.exposure"))

#Checking that X is cool:

head(X)

acc425[which(acc425$chr_pos == "1:693731"),]
WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$chr_pos == "1:693731"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

WholeBodyFatFreeMass_match <- WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$chr_pos%in%X$snp),]

WholeBodyFatFreeMass_match <- WholeBodyFatFreeMass_match[order(match(WholeBodyFatFreeMass_match$chr_pos, X$snp)),]

which(WholeBodyFatFreeMass_match$chr_pos != X$snp) #they all match. 

head(WholeBodyFatFreeMass_match$chr_pos) #perfect
head(X$snp) #perfect

#We have 530 variants that do not possess rsIDs in Body Fat % GWAS.

check <- WholeBodyFatFreeMass_match[order(WholeBodyFatFreeMass_match$SNP),]

head(check, 300) #We do not have the rsid problem here!!

#Let's check if that is true...

index_clean <- which(str_detect(WholeBodyFatFreeMass_match$SNP, "rs") != FALSE)

X_clean <- X[index_clean,]
WholeBodyFatFreeMass_clean <- WholeBodyFatFreeMass_match[index_clean,]

which(X_clean$snp != WholeBodyFatFreeMass_clean$chr_pos) #perfect

head(X_clean$snp) #perfect
head(WholeBodyFatFreeMass_clean$chr_pos) #perfect

tail(X_clean$snp)
tail(WholeBodyFatFreeMass_clean$chr_pos)

X_clean$snp <- WholeBodyFatFreeMass_clean$SNP #now we can return the RSID.

X_end <- X_clean

head(X_end) #perfect
tail(X_end) #perfect

#Calculating Nuisance:

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_acc425_WholeBodyFatFreeMass_ORIGINAL.rds")
saveRDS(params, file = "~/params_acc425_WholeBodyFatFreeMass_ORIGINAL.rds")

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

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/LDShrink/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

WholeBodyFatFreeMass_pruned <- WholeBodyFatFreeMass_match[which(WholeBodyFatFreeMass_match$SNP%in%pruned_all),]

acc425_pruned <- acc425[which(acc425$chr_pos%in%WholeBodyFatFreeMass_pruned$chr_pos),]

write.table(acc425_pruned, "~/acc425_pruned_2_WholeBodyFatFreeMass.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Let's retrieve the SNPs:

top_acc425_pruned_df <- fread("~/acc425_pruned_2_WholeBodyFatFreeMass.txt")

top_acc425_pruned_snps <- top_acc425_pruned_df$SNP

#Let's calcualte this!

res <- cause(X=X_end, variants = top_acc425_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_acc425_to_WholeBodyFatFreeMass.rds")

res <- readRDS("~/RES_acc425_to_WholeBodyFatFreeMass.rds")

#There seem to be problems. Let's check:

res$elpd

#model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing  -3.528542     1.9510749 -1.808512
#2    null  causal  -4.834766     2.7382631 -1.765632
#3 sharing  causal  -1.306224     0.7912341 -1.650870

summary(res)

#p-value testing that causal model is a better fit:  0.049 
#Posterior medians and  95 % credible intervals:
#  model     gamma             eta                    q                  
#[1,] "Sharing" NA                "-0.07 (-0.16, -0.02)" "0.58 (0.08, 0.93)"
#[2,] "Causal"  "-0.05 (-0.1, 0)" "0 (-0.25, 0.21)"      "0.18 (0, 0.85)"   

plot(res, type="data")