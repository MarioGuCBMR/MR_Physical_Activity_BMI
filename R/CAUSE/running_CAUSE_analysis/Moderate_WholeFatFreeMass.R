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
Moderate<- fread("~/RAW_DATA/Doherty-2018-NatureComms-moderate.csv.gz")

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
WholeBodyFatFreeMass$chr_pos <- paste(WholeBodyFatFreeMass$chr.exposure, WholeBodyFatFreeMass$pos.exposure, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(WholeBodyFatFreeMass, Moderate, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("beta.exposure", "BETA"), 
                se_cols = c("se.exposure", "SE"), 
                A1_cols = c("effect_allele.exposure", "ALLELE1"), 
                A2_cols = c("other_allele.exposure", "ALLELE0"))


#Checking that X is cool:

head(X)

Moderate[which(Moderate$chr_pos == "1:102884223"),]
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
#our dear data:

Moderate_not_clean <- Moderate[which(Moderate$chr_pos%in%X_not_clean$snp),]

Moderate_not_clean <- Moderate_not_clean[order(match(Moderate_not_clean$chr_pos, X_not_clean$snp)),]

#Let's check...

head(Moderate_not_clean)
tail(Moderate_not_clean) #some of them can be recovered!!!

index_Moderate_clean <- which(str_detect(Moderate_not_clean$SNP, "rs") == TRUE) #none of them can be recovered with ST....

##################################################
#Let's try to retrieve the reModeratewith PhenoScanner#
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

X_not_clean$snp <- Moderate_not_clean$SNP #they still have the same chr:pos, so the weird ones will be replaced by the rsids

X_end <- rbind(X_clean, X_not_clean)

head(X_end)
tail(X_end)

#Calculating Nuisance:

set.seed(100)
varlist<- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_WholeBodyFatFreeMass_Moderate.rds")
saveRDS(params, file = "~/params_WholeBodyFatFreeMass_Moderate.rds")

#Params has been calculated, so let's go for LD pruning if possible:

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

pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

WholeBodyFatFreeMass_pruned <- WholeBodyFatFreeMass_match[which(WholeBodyFatFreeMass_match$SNP%in%pruned_all),]

write.table(WholeBodyFatFreeMass_pruned, "~/WholeBodyFatFreeMass_pruned_Moderate.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

######################################################

top_WholeBodyFatFreeMass_pruned_df <- read.table("~/WholeBodyFatFreeMass_pruned_Moderate.txt", header = TRUE, stringsAsFactors = FALSE)
top_WholeBodyFatFreeMass_pruned_snps <- top_WholeBodyFatFreeMass_pruned_df$SNP

#Let's run it!!

res <- cause(X=X_end, variants = top_WholeBodyFatFreeMass_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_WholeBodyFatFreeMass_to_Moderate.rds")

res <- readRDS("~/RES_WholeBodyFatFreeMass_to_Moderate.rds")

res$elpd

#model1  model2  delta_elpd se_delta_elpd          z
#1    null sharing -80.7742104    12.7688397 -6.3258849
#2    null  causal -80.8941720    12.9467071 -6.2482430
#3 sharing  causal  -0.1199617     0.7405333 -0.1619936

summary(res)

#p-value testing that causal model is a better fit:  0.44 
#Posterior medians and  95 % credible intervals:
#  model     gamma              eta                    q                  
#[1,] "Sharing" NA                 "-0.19 (-0.32, -0.12)" "0.59 (0.35, 0.89)"
#[2,] "Causal"  "-0.11 (-0.28, 0)" "-0.04 (-0.5, 0.36)"   "0.22 (0.01, 0.85)"

plot(res, type="data")

################################################################################
################################################################################
################################################################################

#################Moderate WholeBodyFatFreeMass################

#########################
#Making space for memory#
#########################

memory.limit(size=800000)

##############
#Loading data#
##############

WholeBodyFatFreeMass <- fread("~/CURATED_DATA/whole_fat_free_mass_UKBB_Combined_Curated_FULL.txt")
Moderate<- fread("~/RAW_DATA/Doherty-2018-NatureComms-moderate.csv.gz")

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
WholeBodyFatFreeMass$chr_pos <- paste(WholeBodyFatFreeMass$chr.exposure, WholeBodyFatFreeMass$pos.exposure, sep = ":")

#######
#CAUSE#
#######

library(cause)

X <- gwas_merge(Moderate, WholeBodyFatFreeMass, snp_name_cols = c("chr_pos", "chr_pos"), 
                beta_hat_cols = c("BETA", "beta.exposure"), 
                se_cols = c("SE", "se.exposure"), 
                A1_cols = c("ALLELE1", "effect_allele.exposure"), 
                A2_cols = c("ALLELE0", "other_allele.exposure"))

#Checking that X is cool:

head(X)

Moderate[which(Moderate$chr_pos == "1:693731"),]
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

Moderate_not_clean <- Moderate[which(Moderate$chr_pos%in%X_not_clean$snp),]

Moderate_not_clean <- Moderate_not_clean[order(match(Moderate_not_clean$chr_pos, X_not_clean$snp)),]

#Let's check...

head(Moderate_not_clean)
tail(Moderate_not_clean) 

index_Moderate_clean <- which(str_detect(Moderate_not_clean$SNP, "rs") == TRUE) #none of them can be recovered with ST....

##################################################
#Let's try to retrieve the reModeratewith PhenoScanner#
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

X_not_clean$snp <- Moderate_not_clean$SNP #they still have the same chr:pos, so the weird ones will be replaced by the rsids

X_end <- rbind(X_clean, X_not_clean)

head(X_end)
tail(X_end)

#Calculating Nuisance:

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "~/X_Moderate_WholeBodyFatFreeMass.rds")
saveRDS(params, file = "~/params_Moderate_WholeBodyFatFreeMass.rds")

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

WholeBodyFatFreeMass_pruned <- WholeBodyFatFreeMass_match[which(WholeBodyFatFreeMass_match$SNP%in%pruned_all),]

Moderate_pruned <- Moderate[which(Moderate$chr_pos%in%WholeBodyFatFreeMass_pruned$chr_pos),]

write.table(Moderate_pruned, "~/Moderate_pruned_2_WholeBodyFatFreeMass_ORIGINAL.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

top_Moderate_pruned_df <- fread("~/Moderate_pruned_2_WholeBodyFatFreeMass_ORIGINAL.txt")

#Let's get the SNPs:
#Careful: in Modeate the variant is still like 9:13982834_G_A, but in WholeBodyFatFreeMass we have
#the variant with rsid. Thus, we are gonna use the rsid from WholeBodyFatFreeMass...

top_Moderate_pruned_snps <- WholeBodyFatFreeMass_match$SNP[which(WholeBodyFatFreeMass_match$chr_pos%in%top_Moderate_pruned_df$chr_pos)]

#Let's calcualte this!

res <- cause(X=X_end, variants = top_Moderate_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "~/RES_Moderate_to_WholeBodyFatFreeMass.rds")

res <- readRDS("~/RES_Moderate_to_WholeBodyFatFreeMass.rds")

#There seem to be problems. Let's check:

res$elpd

#model1  model2 delta_elpd se_delta_elpd         z
#1    null sharing  -3.723475     2.1815940 -1.706768
#2    null  causal  -4.843199     2.8796423 -1.681875
#3 sharing  causal  -1.119725     0.7065897 -1.584689

summary(res)

#p-value testing that causal model is a better fit:  0.057 
#Posterior medians and  95 % credible intervals:
#model     gamma              eta                   q                  
#[1,] "Sharing" NA                 "-0.1 (-0.24, -0.03)" "0.56 (0.09, 0.92)"
#[2,] "Causal"  "-0.07 (-0.15, 0)" "0 (-0.34, 0.31)"     "0.19 (0, 0.86)" 

plot(res, type="data")