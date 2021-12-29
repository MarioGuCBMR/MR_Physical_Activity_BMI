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

WholeBodyFatFreeMass <- fread("~/CURATED_DATA/whole_fat_free_mass_UKBB_Combined_Curated_FULL.txt")
ST <- fread("~/RAW_DATA/Doherty-2018-NatureComms-sedentary.csv.gz")

###########################
#Cleaning WholeBodyFatFreeMass data#
###########################

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
WholeBodyFatFreeMass$chr_pos <- paste(WholeBodyFatFreeMass$chr.exposure, WholeBodyFatFreeMass$pos.exposure, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(WholeBodyFatFreeMass, ST, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("beta.exposure", "BETA"), 
                se_cols = c("se.exposure", "SE"), 
                A1_cols = c("effect_allele.exposure", "ALLELE1"), 
                A2_cols = c("other_allele.exposure", "ALLELE0"))


#Checking that X is cool:

head(X)

ST[which(ST$chr_pos == "1:102884223"),]
WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$chr_pos == "1:102884223"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

WholeBodyFatFreeMass_match <- WholeBodyFatFreeMass[which(WholeBodyFatFreeMass$chr_pos%in%X$snp),]

WholeBodyFatFreeMass_match <- WholeBodyFatFreeMass_match[order(match(WholeBodyFatFreeMass_match$chr_pos, X$snp)),]

which(WholeBodyFatFreeMass_match$chr_pos != X$snp) #they all match. 

#We have 530 variants that do not possess rsIDs.

check <- WholeBodyFatFreeMass_match[order(WholeBodyFatFreeMass_match$SNP),]

head(check, 300) #there are 300~ variants without RSID. Let's try to obtain them...

#######################################
#First we are gonna get the clean ones#
#######################################

index_clean <- which(str_detect(WholeBodyFatFreeMass_match$SNP, "rs") != FALSE)

X_clean <- X[index_clean,]
WholeBodyFatFreeMass_clean <- WholeBodyFatFreeMass_match[index_clean,]

which(X_clean$snp != WholeBodyFatFreeMass_clean$chr_pos)

head(X_clean$snp)
head(WholeBodyFatFreeMass_clean$chr_pos)

tail(X_clean$snp)
tail(WholeBodyFatFreeMass_clean$chr_pos)

X_clean$snp <- WholeBodyFatFreeMass_clean$SNP #now we can return the RSID.

##########
#All good#
##########

index_not_clean <- which(str_detect(WholeBodyFatFreeMass_match$SNP, "rs") == FALSE)

X_not_clean <- X[index_not_clean,]

#And now we are gonna retrieve the chromosome and position from 
#our dear STdata:

ST_not_clean <- ST[which(ST$chr_pos%in%X_not_clean$snp),]

ST_not_clean <- ST_not_clean[order(match(ST_not_clean$chr_pos, X_not_clean$snp)),]

#Let's check...

head(ST_not_clean)
tail(ST_not_clean) #some of them can be recovered!!!

index_ST_clean <- which(str_detect(ST_not_clean$SNP, "rs") == TRUE) #none of them can be recovered with ST....

##################################################
#Let's try to retrieve the rest with PhenoScanner#
##################################################

X_not_clean$snp_ <- paste("chr", X_not_clean$snp, sep = "")

#Finally, let's check in PhenoScanner if these exist:

check_ps_1 <- phenoscanner::phenoscanner(X_not_clean$snp_[seq(1,100)])
check_ps_2 <- phenoscanner::phenoscanner(X_not_clean$snp_[seq(101,200)])
check_ps_3 <- phenoscanner::phenoscanner(X_not_clean$snp_[seq(201,257)])

results <- rbind(check_ps_1$snps, check_ps_2$snps, check_ps_3$snps) #we recovered two! 

X_not_clean[which(X_not_clean$snp_%in%results$snp),] #the alleles do not match!! In PS they are insertions/deletions

################################################
#I checked in dbSNP and the SNPs do exist there#
################################################

#We will have to make do with this. 
#We will still introduce the non-clean, with their original IDS.
#Maybe CAUSE presents these IDs too, you never now.

X_not_clean <- X_not_clean %>%
  select(-snp_)

X_not_clean$snp <- ST_not_clean$SNP #they still have the same chr:pos, so the weird ones will be replaced by the rsids

X_end <- rbind(X_clean, X_not_clean)

head(X_end)
tail(X_end)

#Calculating Nuisance:

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_WholeBodyFatFreeMass_ST.rds")
saveRDS(params, file = "~/params_WholeBodyFatFreeMass_ST.rds")

#ST params has been calculated, so let's go for LD pruning if possible:

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

WholeBodyFatFreeMass_pruned <- WholeBodyFatFreeMass_match[which(WholeBodyFatFreeMass_match$SNP%in%pruned_all),]

write.table(WholeBodyFatFreeMass_pruned, "~/WholeBodyFatFreeMass_ST_pruned.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

######################################################

top_WholeBodyFatFreeMass_pruned_df <- read.table("~/WholeBodyFatFreeMass_ST_pruned.txt", header = TRUE, stringsAsFactors = FALSE)
top_WholeBodyFatFreeMass_pruned_snps <- top_WholeBodyFatFreeMass_pruned_df$SNP

#Let's run it!!

res <- cause(X=X_end, variants = top_WholeBodyFatFreeMass_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_WholeBodyFatFreeMass_to_ST.rds")

res <- readRDS("~/RES_WholeBodyFatFreeMass_to_ST.rds")

res$elpd

#model1  model2  delta_elpd se_delta_elpd         z
#1    null sharing -96.5233253    13.3403448 -7.235445
#2    null  causal -97.4162918    13.4826160 -7.225326
#3 sharing  causal  -0.8929665     0.3619247 -2.467272

summary(res)

#p-value testing that causal model is a better fit:  0.0068 
#Posterior medians and  95 % credible intervals:
#  model     gamma              eta                 q                  
#[1,] "Sharing" NA                 "0.16 (0.13, 0.21)" "0.86 (0.65, 0.98)"
#[2,] "Causal"  "0.14 (0.08, 0.2)" "0 (-0.27, 0.36)"   "0.18 (0, 0.83)"   

plot(res, type="data")

################################################################################
################################################################################
################################################################################

#################ST to WholeBodyFatFreeMass################

#########################
#Making space for memory#
#########################

memory.limit(size=800000)

##############
#Loading data#
##############

WholeBodyFatFreeMass <- fread("~/CURATED_DATA/whole_fat_free_mass_UKBB_Combined_Curated_FULL.txt")
ST <- fread("~/RAW_DATA/Doherty-2018-NatureComms-sedentary.csv.gz")

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
WholeBodyFatFreeMass$chr_pos <- paste(WholeBodyFatFreeMass$chr.exposure, WholeBodyFatFreeMass$pos.exposure, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(ST, WholeBodyFatFreeMass, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "beta.exposure"), 
                se_cols = c("SE", "se.exposure"), 
                A1_cols = c("ALLELE1", "effect_allele.exposure"), 
                A2_cols = c("ALLELE0", "other_allele.exposure"))

#Checking that X is cool:

head(X)

ST[which(ST$chr_pos == "1:693731"),]
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

head(check, 300) #there are 300~ variants without RSID. Let's try to obtain them...

#######################################
#First we are gonna get the clean ones#
#######################################

index_clean <- which(str_detect(WholeBodyFatFreeMass_match$SNP, "rs") != FALSE)

X_clean <- X[index_clean,]
WholeBodyFatFreeMass_clean <- WholeBodyFatFreeMass_match[index_clean,]

which(X_clean$snp != WholeBodyFatFreeMass_clean$chr_pos) #perfect

head(X_clean$snp) #perfect
head(WholeBodyFatFreeMass_clean$chr_pos) #perfect

tail(X_clean$snp)
tail(WholeBodyFatFreeMass_clean$chr_pos)

X_clean$snp <- WholeBodyFatFreeMass_clean$SNP #now we can return the RSID.

##########
#All good#
##########

index_not_clean <- which(str_detect(WholeBodyFatFreeMass_match$SNP, "rs") == FALSE)

X_not_clean <- X[index_not_clean,]

#And now we are gonna retrieve the chromosome and position from 
#our dear STdata:

ST_not_clean <- ST[which(ST$chr_pos%in%X_not_clean$snp),]

ST_not_clean <- ST_not_clean[order(match(ST_not_clean$chr_pos, X_not_clean$snp)),]

#Let's check...

head(ST_not_clean)
tail(ST_not_clean) 

index_ST_clean <- which(str_detect(ST_not_clean$SNP, "rs") == TRUE) #none of them can be recovered with ST....

##################################################
#Let's try to retrieve the rest with PhenoScanner#
##################################################

X_not_clean$snp_ <- paste("chr", X_not_clean$snp, sep = "")

#Finally, let's check in PhenoScanner if these exist:

check_ps_1 <- phenoscanner::phenoscanner(X_not_clean$snp_[seq(1,100)])
check_ps_2 <- phenoscanner::phenoscanner(X_not_clean$snp_[seq(101,200)])
check_ps_3 <- phenoscanner::phenoscanner(X_not_clean$snp_[seq(201,257)])

results <- rbind(check_ps_1$snps, check_ps_2$snps, check_ps_3$snps) #we recovered two! 

X_not_clean[which(X_not_clean$snp_%in%results$snp),] #the alleles do not match!! In PS they are insertions/deletions

################################################
#I checked in dbSNP and the SNPs do exist there#
################################################

#We will have to make do with this. 
#We will still introduce the non-clean, with their original IDS.
#Maybe CAUSE presents these IDs too, you never now.

X_not_clean <- X_not_clean %>%
  select(-snp_)

X_not_clean$snp <- ST_not_clean$SNP #they still have the same chr:pos, so the weird ones will be replaced by the rsids

X_end <- rbind(X_clean, X_not_clean)

head(X_end)
tail(X_end)

#Calculating Nuisance:

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_ST_WholeBodyFatFreeMass.rds")
saveRDS(params, file = "~/params_ST_WholeBodyFatFreeMass.rds")

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

WholeBodyFatFreeMass_pruned <- WholeBodyFatFreeMass_match[which(WholeBodyFatFreeMass_match$SNP%in%pruned_all),]

ST_pruned <- ST[which(ST$chr_pos%in%WholeBodyFatFreeMass_pruned$chr_pos),]

write.table(ST_pruned, "~/ST_WholeBodyFatFreeMass_pruned.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

top_ST_pruned_df <- fread("~/ST_WholeBodyFatFreeMass_pruned.txt")

#Let's get the SNPs:

top_ST_pruned_snps <- top_ST_pruned_df$SNP

#Let's calcualte this!

res <- cause(X=X_end, variants = top_ST_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_ST_to_WholeBodyFatFreeMass.rds")

res <- readRDS("~/RES_ST_to_WholeBodyFatFreeMass.rds")

#There seem to be problems. Let's check:

res$elpd

#model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing  -3.099351     1.7911589 -1.730361
#2    null  causal  -4.295064     2.5122419 -1.709654
#3 sharing  causal  -1.195713     0.7231511 -1.653476

summary(res)

#p-value testing that causal model is a better fit:  0.049 
#Posterior medians and  95 % credible intervals:
#  model     gamma            eta                 q                  
#[1,] "Sharing" NA               "0.07 (0.02, 0.17)" "0.56 (0.06, 0.92)"
#[2,] "Causal"  "0.06 (0, 0.11)" "0 (-0.23, 0.23)"   "0.19 (0, 0.85)"   

plot(res, type="data")