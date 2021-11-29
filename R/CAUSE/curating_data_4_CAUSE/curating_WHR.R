##############
#INTRODUCTION#
##############

#This is a quick code to check how much do we need to clean the whr stratified data.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)
library(TwoSampleMR)

##################
#Loading function#
##################

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

clean_liftover <- function(chr_pos){
  
  tmp_1 <- strsplit(chr_pos, ":")[[1]][1]
  tmp_2 <- strsplit(chr_pos, ":")[[1]][2]
  tmp_3 <- as.numeric(as.character(strsplit(tmp_2, "-")[[1]][1])) +1
  
  final <- paste(tmp_1, tmp_3, sep = ":")
  
  return(final)
  
}

##############
#Loading data#
##############

whr <- fread("~/RAW_DATA/GIANT_2015_WHR_COMBINED_EUR.txt") 

#Let's check what columns do we have:

whr <- whr[order(whr$MarkerName),]

head(whr, 1000)

#Here we just have the chromsomes and the positions...
#But we do not have the MAF so we will remove them.

##############################
#Filtering data and saving it#
##############################

#What we are going to do is the following:

#a) Filter for MAF/INFO if possible.
#b) remove MHC region SNPs here, if possible. If not, in the outcome.
#c) Save the data so we do not have to wait for the long query in the database.

whr_maf <- whr[which(whr$FreqAllele1HapMapCEU > 0.01),]
whr_maf <- whr_maf[which(whr_maf$FreqAllele1HapMapCEU < 0.99),]

summary(whr_maf$FreqAllele1HapMapCEU) #perfect.

#Now the INFO, do we have that data?

colnames(whr_maf) #we do not. 

##########################
#REMOVING MHC region SNPs#
##########################

#WAIT. We cannot do this because we use the filters from build 37 from Warrington et al.
#Hence, what we are gonna do is... just include them.
#Since we removed them in the exposure, they will be removed automatically.

################################################
#Getting as many chr_pos possible from BMI data#
################################################

BMI <- fread("~/RAW_DATA/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")

BMI$RSID <- as.character(unlist(sapply(BMI$SNP, parse_rsid)))

BMI <- BMI[which(is.na(BMI$INFO) == FALSE),]
BMI <- BMI[which(is.na(BMI$CHR) == FALSE),]

summary(BMI$CHR) #perfect
summary(BMI$INFO) #perfect

#We removed the weird SNPs.

chr_ <- paste("chr", BMI$CHR, sep = "")
BMI$chr_pos <- paste(chr_, BMI$POS, sep = ":")

#Now we are gonna match by SNP, 
#but we might miss some of them because they were merged.
#But we can do that later.

BMI_match <- BMI[which(BMI$RSID%in%whr_maf$MarkerName),] #most of them:

whr_Pulit <- whr_maf[which(whr_maf$MarkerName%in%BMI_match$RSID),]

#We have duplicates, but they are bound to have the same chr_pos, so we do not care
#about them. Let's check just in case:

dupl <- BMI_match$RSID[which(duplicated(BMI_match$RSID) == TRUE)]

BMI_dupl <- BMI_match[which(BMI_match$RSID%in%dupl),]

#Let's order them and get ready the dataframe:

BMI_dupl <- BMI_dupl[order(BMI_dupl$RSID),]

head(BMI_dupl) #same chr_pos
tail(BMI_dupl) #same chr_pos

#PERFECT. We can remove duplicates:

BMI_match <- BMI_match[which(duplicated(BMI_match$RSID) == FALSE),] #now they match.

BMI_match <- BMI_match[order(match(BMI_match$RSID, whr_Pulit$MarkerName)),]

#Let's check if we did it properly:

which(BMI_match$RSID != whr_Pulit$MarkerName) #perfect.

head(BMI_match$RSID)
head(whr_Pulit$MarkerName)

#We can match the chr_pos:

whr_Pulit$chr_pos_37 <- BMI_match$chr_pos

#########
#PERFECT#
#########

#Let's go and get those missing:

whr_missing <- whr_maf[which(!(whr_maf$MarkerName%in%BMI_match$RSID)),]

#Let's see how many of them can be recovered from SNP Nexus:

for_nexus <- whr_missing

for_nexus$dbsnp <- "dbsnp"

for_nexus <- for_nexus %>%
  select(dbsnp, MarkerName)

for_nexus <- for_nexus[order(for_nexus$MarkerName),]

head(for_nexus)
tail(for_nexus)

########################################################
#Important we have three SNP that are in another format#
########################################################

#We are gonna use liftover to transform them.
#Though technically there is no need to,
#They are the same as we had in WHRAdjBMI combined.

#chr1:161981977-161981979
#chr6:13159524-13159526
#chr6:26020935-26020937

gotcha_vect <- c("chr1:161981978", "chr6:13159525", "chr6:26020936")

which(BMI$chr_pos%in%gotcha_vect) #can we recover them with BMI? YES.

BMI_gotcha <- BMI[which(BMI$chr_pos%in%gotcha_vect),]

###############################################
#Let's prepare the data for getting these snps#
###############################################

whr_missing <- whr_missing[order(whr_missing$MarkerName),]

whr_missing_gotcha <- whr_missing[seq(1,3),]

#After checking, the SNP match. 
#Hence:

whr_missing_gotcha$MarkerName <- BMI_gotcha$SNP

whr_missing_gotcha$chr_pos_37 <- BMI_gotcha$chr_pos

#########
#PERFECT#
#########

#With that being settled..., we only need to run SNPNexus:

head(for_nexus)

for_nexus <- for_nexus[-seq(1,3),]

head(for_nexus) #PERFECT:

write.table(for_nexus, "whr_combined_4_Nexus.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

####################
#Let's recover them#
####################

recovered_snps <- fread("whr_combined_in_Nexus.txt")

###################################
#We are gonna be careful and check#
###################################

recovered_snps <- recovered_snps[order(recovered_snps$`REF Allele`),]

head(recovered_snps) #one deletion
tail(recovered_snps) #nothing.

#It seems we only need to remove the deletion. 
#To do so, since we have IUPAC alleles we will just take all the alleles
#that are letters, to avoid insertions and deletions.

good_alleles <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
                  "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
                  "V", "W", "X", "Y", "Z")

recovered_snps <- recovered_snps[which(recovered_snps$`REF Allele`%in%good_alleles),] #19511->19510. PERFECT.
recovered_snps <- recovered_snps[which(recovered_snps$`ALT Allele (IUPAC)`%in%good_alleles),] #19510->19510. PERFECT.

##########
#ALL GOOD#
##########

whr_recovered <- whr_missing[which(whr_missing$MarkerName%in%recovered_snps$dbSNP),]

recovered_snps <- recovered_snps[order(match(recovered_snps$dbSNP, whr_recovered$MarkerName)),]

#Doing some checkity checks:

which(recovered_snps$dbSNP != whr_recovered$MarkerName)

head(recovered_snps$dbSNP)
head(whr_recovered$MarkerName)

#PERFECT:

chr_ <- paste("chr", recovered_snps$Chromosome, sep = "")

recovered_snps$chr_pos <- paste(chr_, recovered_snps$Position, sep = ":")

whr_recovered$chr_pos_37 <- recovered_snps$chr_pos

###############################
#Let's check the missing ones:#
###############################

whr_missing_2 <- whr_missing[which(!(whr_missing$MarkerName%in%recovered_snps$dbSNP)),]

View(whr_missing_2)

#The first three were already recovered, so let's remove them:

whr_missing_2 <- whr_missing_2[-seq(1,3),]

head(whr_missing_2)

#######################################################################
#Let's do one more run of liftover to recover as many SNPs as possible#
#######################################################################

chr_ <- paste("chr", whr_missing_2$Chr, sep = "")

pos <- paste(whr_missing_2$Pos -1, "-", whr_missing_2$Pos+1, sep = "")

chr_pos_18 <- paste(chr_, pos, sep = ":")

whr_missing_2$chr_pos_18 <- chr_pos_18

#############################
#Okay, let's do the liftover#
#############################

write.table(chr_pos_18, "/WHR_4_liftover.txt", quote = FALSE, row.names = FALSE, col.names =  FALSE)

#########################################################################################
#Now we are gonna redo Liftover without the failing conversions so we can match the data#
#########################################################################################

failed_conversions <- fread("Failed_Conversions.txt")

whr_missing_2_new <- whr_missing_2[which(!(whr_missing_2$chr_pos_18%in%failed_conversions$`#Deleted in new`)),]

#The numbers match: we had 175 failed conversions.

write.table(whr_missing_2_new$chr_pos_18, "WHR_4_liftover_curated.txt", quote = FALSE, row.names = FALSE, col.names =  FALSE)

converted_data <- as.data.frame(read.table("whr_liftover_correct.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

##########################
#Let's get the chr_pos_37#
##########################

chr_pos_37 <- as.character(unlist(sapply(converted_data$V1, clean_liftover)))

whr_missing_2_new$chr_pos_37 <- chr_pos_37

head(whr_missing_2_new)
tail(whr_missing_2_new)

#Perfect!!!

######################################################
#Now... what are the ones that we could not retrieve?#
######################################################

whr_missing_3 <- whr_missing_2[which(!(whr_missing_2$MarkerName%in%whr_missing_2_new$MarkerName)),]

#Let's take a quick look at them:

whr_missing_3[order(whr_missing_3$MarkerName),] #they present RSIDs... 

#Let's try with phenoscanner. Our last resort:

results_1 <- phenoscanner::phenoscanner(whr_missing_3$MarkerName[seq(1,100)])
results_2 <- phenoscanner::phenoscanner(whr_missing_3$MarkerName[seq(100,175)])

results <- rbind(results_1$snps, results_2$snps) #0 SNPs found in phenoscanner. Maybe because it is down.

################################
#Are some of the SNPs merged???#
################################

merged_rs <- fread("~/RAW_DATA/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

index_merged <- which(merged_rs$old_rs%in%whr_missing_3$MarkerName) #2 of them:

swapping_snps <- merged_rs[index_merged,]

swapping_snps

#V1        V2  V3 V4                      V5                      V6        V7 V8            V9
#1:      668 281865545 136  1 2014-08-25 17:55:49.607 2014-08-25 17:55:49.607 281865545  1 JIRA SNP-6845
#2: 12020931   9604473 131  0   2009-12-02 15:53:00.0   2010-03-15 13:21:00.0   9604473  0           rsm
#old_rs      new_rs
#1:      rs668 rs281865545
#2: rs12020931   rs9604473

#I checked them manually.
#They all have weird positions in build 37 that do not allow to locate them properly.
#So far, we are just gonna give them the "-" treatment. If the RSIDs are found while matching...
#Then better.

whr_missing_3$chr_pos_37 <- "-"

#############
#CONCLUSIONS#
#############

#We solved most of the SNPs, except for 175. Not bad...
#whr_missing_2 and whr_missing_3 have data for chr_pos_18, the others don't.
#so let's make dummy variables to match the df:

whr_Pulit$chr_pos_18 <- "-"
whr_recovered$chr_pos_18 <- "-"
whr_missing_gotcha$chr_pos_18 <- "-"

whr_end <- rbind(whr_Pulit, whr_recovered, whr_missing_gotcha, whr_missing_2_new, whr_missing_3) #PERFECT number. It matches the data perfectily.

dim(whr_end)
dim(whr_maf)

fwrite(whr_end, "~/CURATED_DATA/WHR_combined_Curated.txt")