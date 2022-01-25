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

BodyFatPerc <- fread("~/CURATED_DATA/body_fat_UKBB_Combined_Curated_FULL.txt")
Moderate<- fread("~/RAW_DATA/Doherty-2018-NatureComms-moderate.csv.gz")

###########################
#Cleaning BodyFatPerc data#
###########################

#(Since we had to transform the data from vcf to a readble txt with all the necessary column)
#all of the changes that we are about to make have already been done!)
#Nonetheless, I will include them either way.

summary(BodyFatPerc$eaf.exposure) 

#We do not have INFO so we skip that..., but we do have alleles! CAUSE will only detect
#standard alleles

yes_vect <- c("A", "G", "C", "T")

BodyFatPerc <- BodyFatPerc[which(BodyFatPerc$effect_allele.exposure%in%yes_vect),]
BodyFatPerc <- BodyFatPerc[which(BodyFatPerc$other_allele.exposure%in%yes_vect),]

#And now we remove the data in the MHC region:

#The expanded MHC is 26M to 24M.

BodyFatPerc_mhc <- BodyFatPerc[which(BodyFatPerc$chr.exposure == 6 & BodyFatPerc$pos.exposure > 26000000 & BodyFatPerc$pos.exposure < 34000000),] #so we are gonna approximate it this way.

#All of them were removed in the curation!!

#Perfect!

#########################
#Cleaning Sedentary data#
#########################

#MAF > 0.01
#INFO > 0.7
#Also, removthing those that are deletions or insertions (classic alleles.)

summary(Moderate$A1FREQ) #maf is already done!
summary(Moderate$INFO) #INFO needs to be done:

Moderate<- Moderate[which(Moderate$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

Moderate<- Moderate[which(Moderate$ALLELE1%in%yes_vect),]
Moderate<- Moderate[which(Moderate$ALLELE0%in%yes_vect),]

########################################
#Matching data in the beModerateway possible#
########################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure firModerateand do the matching later to save some of them.

Moderate$chr_pos <- paste(Moderate$CHR, Moderate$BP, sep = ":")
BodyFatPerc$chr_pos <- paste(BodyFatPerc$chr.exposure, BodyFatPerc$pos.exposure, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(BodyFatPerc, Moderate, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("beta.exposure", "BETA"), 
                se_cols = c("se.exposure", "SE"), 
                A1_cols = c("effect_allele.exposure", "ALLELE1"), 
                A2_cols = c("other_allele.exposure", "ALLELE0"))


#Checking that X is cool:

head(X)

Moderate[which(Moderate$chr_pos == "1:102884223"),]
BodyFatPerc[which(BodyFatPerc$chr_pos == "1:102884223"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

BodyFatPerc_match <- BodyFatPerc[which(BodyFatPerc$chr_pos%in%X$snp),]

BodyFatPerc_match <- BodyFatPerc_match[order(match(BodyFatPerc_match$chr_pos, X$snp)),]

which(BodyFatPerc_match$chr_pos != X$snp) #they all match. 

#We have 530 variants that do not possess rsIDs.

check <- BodyFatPerc_match[order(BodyFatPerc_match$SNP),]

head(check, 300) #there are 300~ variants without RSID. Let's try to obtain them...

#######################################
#First we are gonna get the clean ones#
#######################################

index_clean <- which(str_detect(BodyFatPerc_match$SNP, "rs") != FALSE)

X_clean <- X[index_clean,]
BodyFatPerc_clean <- BodyFatPerc_match[index_clean,]

which(X_clean$snp != BodyFatPerc_clean$chr_pos)

head(X_clean$snp)
head(BodyFatPerc_clean$chr_pos)

tail(X_clean$snp)
tail(BodyFatPerc_clean$chr_pos)

X_clean$snp <- BodyFatPerc_clean$SNP #now we can return the RSID.

##########
#All good#
##########

index_not_clean <- which(str_detect(BodyFatPerc_match$SNP, "rs") == FALSE)

X_not_clean <- X[index_not_clean,]

#And now we are gonna retrieve the chromosome and position from 
#our dear STdata:

Moderate_not_clean <- Moderate[which(Moderate$chr_pos%in%X_not_clean$snp),]

Moderate_not_clean <- Moderate_not_clean[order(match(Moderate_not_clean$chr_pos, X_not_clean$snp)),]

#Let's check...

head(Moderate_not_clean)
tail(Moderate_not_clean) #some of them can be recovered!!!

index_Moderate_clean <- which(str_detect(Moderate_not_clean$SNP, "rs") == TRUE) #none of them can be recovered with ST....

#######################################################
#Let's try to retrieve the reModeratewith PhenoScanner#
#######################################################

X_not_clean$snp_ <- paste("chr", X_not_clean$snp, sep = "")

#Finally, let's check in PhenoScanner if these exist:

check_ps_1 <- phenoscanner::phenoscanner(X_not_clean$snp_[seq(1,100)])
check_ps_2 <- phenoscanner::phenoscanner(X_not_clean$snp_[seq(101,200)])
check_ps_3 <- phenoscanner::phenoscanner(X_not_clean$snp_[seq(201,257)])

results <- rbind(check_ps_1$snps, check_ps_2$snps, check_ps_3$snps) #we recovered two! 

X_not_clean[which(X_not_clean$snp_%in%results$snp),] #the alleles do not match!! In PS they are insertions/deletions

################################################
#I checked in dbSNP and the SNPs do exiModeratethere#
################################################

#We will have to make do with this. 
#We will still introduce the non-clean, with their original IDS.
#Maybe CAUSE presents these IDs too, you never now.

X_not_clean <- X_not_clean %>%
  select(-snp_)

X_not_clean$snp <- Moderate_not_clean$SNP #they still have the same chr:pos, so the weird ones will be replaced by the rsids

X_end <- rbind(X_clean, X_not_clean)

head(X_end)
tail(X_end)

#Calculating Nuisance:

set.seed(100)
varlist<- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_BodyFatPerc_Moderate.rds")
saveRDS(params, file = "~/params_BodyFatPerc_Moderate.rds")

#####################
#Let' clump the data#
#####################

library(tidyverse)

variants <- X_end %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

#1

ld <- readRDS("~/LD_panels/chr1_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr1_AF0.05_snpdata.RDS")

pruned_1 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#2

ld <- readRDS("~/LD_panels/chr2_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr2_AF0.05_snpdata.RDS")

pruned_2 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#3

ld <- readRDS("~/LD_panels/chr3_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr3_AF0.05_snpdata.RDS")

pruned_3 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#4

ld <- readRDS("~/LD_panels/chr4_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr4_AF0.05_snpdata.RDS")

pruned_4 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#5

ld <- readRDS("~/LD_panels/chr5_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr5_AF0.05_snpdata.RDS")

pruned_5 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#6

ld <- readRDS("~/LD_panels/chr6_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr6_AF0.05_snpdata.RDS")

pruned_6 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#7

ld <- readRDS("~/LD_panels/chr7_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr7_AF0.05_snpdata.RDS")

pruned_7 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#8

ld <- readRDS("~/LD_panels/chr8_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr8_AF0.05_snpdata.RDS")

pruned_8 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#9

ld <- readRDS("~/LD_panels/chr9_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr9_AF0.05_snpdata.RDS")

pruned_9 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#10

ld <- readRDS("~/LD_panels/chr10_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr10_AF0.05_snpdata.RDS")

pruned_10 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#11

ld <- readRDS("~/LD_panels/chr11_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr11_AF0.05_snpdata.RDS")

pruned_11 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#12

ld <- readRDS("~/LD_panels/chr12_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr12_AF0.05_snpdata.RDS")

pruned_12 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#13

ld <- readRDS("~/LD_panels/chr13_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr13_AF0.05_snpdata.RDS")

pruned_13 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#14

ld <- readRDS("~/LD_panels/chr14_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr14_AF0.05_snpdata.RDS")

pruned_14 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#15

ld <- readRDS("~/LD_panels/chr15_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr15_AF0.05_snpdata.RDS")

pruned_15 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#16

ld <- readRDS("~/LD_panels/chr16_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr16_AF0.05_snpdata.RDS")

pruned_16 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#17

ld <- readRDS("~/LD_panels/chr17_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr17_AF0.05_snpdata.RDS")

pruned_17 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#18

ld <- readRDS("~/LD_panels/chr18_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr18_AF0.05_snpdata.RDS")

pruned_18 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#19

ld <- readRDS("~/LD_panels/chr19_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr19_AF0.05_snpdata.RDS")

pruned_19 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#20

ld <- readRDS("~/LD_panels/chr20_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr20_AF0.05_snpdata.RDS")

pruned_20 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#21

ld <- readRDS("~/LD_panels/chr21_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr21_AF0.05_snpdata.RDS")

pruned_21 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))
#22

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "~/LD_panels/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="~/LD_panels/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

BodyFatPerc_pruned <- BodyFatPerc_match[which(BodyFatPerc_match$SNP%in%pruned_all),]

write.table(BodyFatPerc_pruned, "~/BodyFatPerc_pruned_Moderate.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

top_BodyFatPerc_pruned_df <- read.table("~/BodyFatPerc_pruned_Moderate.txt", header = TRUE, stringsAsFactors = FALSE)
top_BodyFatPerc_pruned_snps <- top_BodyFatPerc_pruned_df$SNP

#Let's run it!!

res <- cause(X=X_end, variants = top_BodyFatPerc_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_BodyFatPerc_to_Moderate.rds")

res <- readRDS("~/RES_BodyFatPerc_to_Moderate.rds")

res$elpd

#model1  model2  delta_elpd se_delta_elpd         z
#1    null sharing -170.462808    17.4595899 -9.763277
#2    null  causal -171.495090    17.6282254 -9.728438
#3 sharing  causal   -1.032282     0.5432155 -1.900317

summary(res)

#p-value testing that causal model is a better fit:  0.029 
#Posterior medians and  95 % credible intervals:
#  model     gamma                  eta                    q                  
#[1,] "Sharing" NA                     "-0.21 (-0.25, -0.18)" "0.91 (0.77, 0.99)"
#[2,] "Causal"  "-0.19 (-0.26, -0.12)" "0 (-0.3, 0.31)"       "0.19 (0, 0.86)"   

plot(res, type="data")

################################################################################
################################################################################
################################################################################

#################Moderate BodyFatPerc################

#########################
#Making space for memory#
#########################

memory.limit(size=800000)

##############
#Loading data#
##############

BodyFatPerc <- fread("~CURATED_DATA/body_fat_UKBB_Combined_Curated_FULL.txt")
Moderate<- fread("~/RAW_DATA/Doherty-2018-NatureComms-moderate.csv.gz")

###########################
#Cleaning BodyFatPerc data#
###########################

#(Since we had to transform the data from vcf to a readble txt with all the necessary column)
#all of the changes that we are about to make have already been done!)
#Nonetheless, I will include them either way.

summary(BodyFatPerc$eaf.exposure) 

#We do not have INFO so we skip that..., but we do have alleles! CAUSE will only detect
#standard alleles

yes_vect <- c("A", "G", "C", "T")

BodyFatPerc <- BodyFatPerc[which(BodyFatPerc$effect_allele.exposure%in%yes_vect),]
BodyFatPerc <- BodyFatPerc[which(BodyFatPerc$other_allele.exposure%in%yes_vect),]

#And now we remove the data in the MHC region:

#The expanded MHC is 26M to 24M.

BodyFatPerc_mhc <- BodyFatPerc[which(BodyFatPerc$chr.exposure == 6 & BodyFatPerc$pos.exposure > 26000000 & BodyFatPerc$pos.exposure < 34000000),] #so we are gonna approximate it this way.

#All of them were removed in the curation!!

#Perfect!

#########################
#Cleaning Sedentary data#
#########################

#MAF > 0.01
#INFO > 0.7
#Also, removthing those that are deletions or insertions (classic alleles.)

summary(Moderate$A1FREQ) #maf is already done!
summary(Moderate$INFO) #INFO needs to be done:

Moderate<- Moderate[which(Moderate$INFO > 0.7),]

yes_vect <- c("A", "G", "C", "T")

Moderate<- Moderate[which(Moderate$ALLELE1%in%yes_vect),]
Moderate<- Moderate[which(Moderate$ALLELE0%in%yes_vect),]

########################################
#Matching data in the best way possible#
########################################

#The matching could be done in many ways, but it is too complicated
#CAUSE solves this by removing all duplicates.
#I guess that a trick to improve a bit is removing the duplicates that present a pvalue lower than the other duplicate.
#in the exposure first and do the matching later to save some of them.

Moderate$chr_pos <- paste(Moderate$CHR, Moderate$BP, sep = ":")
BodyFatPerc$chr_pos <- paste(BodyFatPerc$chr.exposure, BodyFatPerc$pos.exposure, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(Moderate, BodyFatPerc, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "beta.exposure"), 
                se_cols = c("SE", "se.exposure"), 
                A1_cols = c("ALLELE1", "effect_allele.exposure"), 
                A2_cols = c("ALLELE0", "other_allele.exposure"))

#Checking that X is cool:

head(X)

Moderate[which(Moderate$chr_pos == "1:693731"),]
BodyFatPerc[which(BodyFatPerc$chr_pos == "1:693731"),]

#the matching of the alleles was perfect.

#Now we are going to change the X snp to rsid so we can 
#filter the heck out of it later.

BodyFatPerc_match <- BodyFatPerc[which(BodyFatPerc$chr_pos%in%X$snp),]

BodyFatPerc_match <- BodyFatPerc_match[order(match(BodyFatPerc_match$chr_pos, X$snp)),]

which(BodyFatPerc_match$chr_pos != X$snp) #they all match. 

head(BodyFatPerc_match$chr_pos) #perfect
head(X$snp) #perfect

#We have 530 variants that do not possess rsIDs in Body Fat % GWAS.

check <- BodyFatPerc_match[order(BodyFatPerc_match$SNP),]

head(check, 300) #there are 300~ variants without RSID. Let's try to obtain them...

#######################################
#First we are gonna get the clean ones#
#######################################

index_clean <- which(str_detect(BodyFatPerc_match$SNP, "rs") != FALSE)

X_clean <- X[index_clean,]
BodyFatPerc_clean <- BodyFatPerc_match[index_clean,]

which(X_clean$snp != BodyFatPerc_clean$chr_pos) #perfect

head(X_clean$snp) #perfect
head(BodyFatPerc_clean$chr_pos) #perfect

tail(X_clean$snp)
tail(BodyFatPerc_clean$chr_pos)

X_clean$snp <- BodyFatPerc_clean$SNP #now we can return the RSID.

##########
#All good#
##########

index_not_clean <- which(str_detect(BodyFatPerc_match$SNP, "rs") == FALSE)

X_not_clean <- X[index_not_clean,]

#And now we are gonna retrieve the chromosome and position from 
#our dear STdata:

Moderate_not_clean <- Moderate[which(Moderate$chr_pos%in%X_not_clean$snp),]

Moderate_not_clean <- Moderate_not_clean[order(match(Moderate_not_clean$chr_pos, X_not_clean$snp)),]

#Let's check...

head(Moderate_not_clean)
tail(Moderate_not_clean) 

index_Moderate_clean <- which(str_detect(Moderate_not_clean$SNP, "rs") == TRUE) #none of them can be recovered with ST....

#######################################################
#Let's try to retrieve the reModeratewith PhenoScanner#
#######################################################

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

X_not_clean$snp <- Moderate_not_clean$SNP #they still have the same chr:pos, so the weird ones will be replaced by the rsids

X_end <- rbind(X_clean, X_not_clean)

head(X_end)
tail(X_end)

#Calculating Nuisance:

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_Moderate_BodyFatPerc.rds")
saveRDS(params, file = "~/params_Moderate_BodyFatPerc.rds")

#Params has been calculated, so let's go for LD pruning if possible:

variants <- X_end %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

#1

ld <- readRDS("~/LD_panels/chr1_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr1_AF0.05_snpdata.RDS")

pruned_1 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#2

ld <- readRDS("~/LD_panels/chr2_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr2_AF0.05_snpdata.RDS")

pruned_2 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#3

ld <- readRDS("~/LD_panels/chr3_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr3_AF0.05_snpdata.RDS")

pruned_3 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#4

ld <- readRDS("~/LD_panels/chr4_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr4_AF0.05_snpdata.RDS")

pruned_4 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#5

ld <- readRDS("~/LD_panels/chr5_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr5_AF0.05_snpdata.RDS")

pruned_5 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#6

ld <- readRDS("~/LD_panels/chr6_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr6_AF0.05_snpdata.RDS")

pruned_6 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#7

ld <- readRDS("~/LD_panels/chr7_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr7_AF0.05_snpdata.RDS")

pruned_7 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#8

ld <- readRDS("~/LD_panels/chr8_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr8_AF0.05_snpdata.RDS")

pruned_8 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#9

ld <- readRDS("~/LD_panels/chr9_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr9_AF0.05_snpdata.RDS")

pruned_9 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#10

ld <- readRDS("~/LD_panels/chr10_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr10_AF0.05_snpdata.RDS")

pruned_10 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#11

ld <- readRDS("~/LD_panels/chr11_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr11_AF0.05_snpdata.RDS")

pruned_11 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#12

ld <- readRDS("~/LD_panels/chr12_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr12_AF0.05_snpdata.RDS")

pruned_12 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#13

ld <- readRDS("~/LD_panels/chr13_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr13_AF0.05_snpdata.RDS")

pruned_13 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#14

ld <- readRDS("~/LD_panels/chr14_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr14_AF0.05_snpdata.RDS")

pruned_14 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#15

ld <- readRDS("~/LD_panels/chr15_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr15_AF0.05_snpdata.RDS")

pruned_15 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#16

ld <- readRDS("~/LD_panels/chr16_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr16_AF0.05_snpdata.RDS")

pruned_16 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#17

ld <- readRDS("~/LD_panels/chr17_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr17_AF0.05_snpdata.RDS")

pruned_17 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#18

ld <- readRDS("~/LD_panels/chr18_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr18_AF0.05_snpdata.RDS")

pruned_18 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#19

ld <- readRDS("~/LD_panels/chr19_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr19_AF0.05_snpdata.RDS")

pruned_19 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#20

ld <- readRDS("~/LD_panels/chr20_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr20_AF0.05_snpdata.RDS")

pruned_20 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#21

ld <- readRDS("~/LD_panels/chr21_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr21_AF0.05_snpdata.RDS")

pruned_21 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))
#22

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "~/LD_panels/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="~/LD_panels/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("~/LD_panels/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("~/LD_panels/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

BodyFatPerc_pruned <- BodyFatPerc_match[which(BodyFatPerc_match$SNP%in%pruned_all),]

Moderate_pruned <- Moderate[which(Moderate$chr_pos%in%BodyFatPerc_pruned$chr_pos),]

write.table(Moderate_pruned, "~/Moderate_pruned_2_BodyFatPerc.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

top_Moderate_pruned_df <- fread("~/Moderate_pruned_2_BodyFatPerc.txt")

#Let's get the SNPs:
#Careful: in Modeate the variant is still like 9:13982834_G_A, but in BodyFatPerc we have
#the variant with rsid. Thus, we are gonna use the rsid from BodyFatPerc...

top_Moderate_pruned_snps <- BodyFatPerc_match$SNP[which(BodyFatPerc_match$chr_pos%in%top_Moderate_pruned_df$chr_pos)]

#Let's calcualte this!

res <- cause(X=X_end, variants = top_Moderate_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_Moderate_to_BodyFatPerc.rds")

res <- readRDS("~/RES_Moderate_to_BodyFatPerc.rds")

#There seem to be problems. Let's check:

res$elpd

#model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing  -8.706653     2.6045635 -3.342845
#2    null  causal -10.108223     3.1324111 -3.226978
#3 sharing  causal  -1.401569     0.5296923 -2.646006

summary(res)

#p-value testing that causal model is a better fit:  0.0041 
#Posterior medians and  95 % credible intervals:
#  model     gamma                  eta                    q                  
#[1,] "Sharing" NA                     "-0.17 (-0.26, -0.09)" "0.78 (0.39, 0.97)"
#[2,] "Causal"  "-0.15 (-0.24, -0.05)" "0 (-0.32, 0.31)"      "0.2 (0, 0.86)" 

plot(res, type="data")