##############
#INTRODUCTION#
##############

#This code is to curated Ever-smoker data from the downloaded vcf file.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(TwoSampleMR)

###################
#Loading functions#
###################

beta_retriever <- function(header_info){
  
  beta <- strsplit(header_info, ":")[[1]][1]
  
  return(beta)
  
}

se_retriever <- function(header_info){
  
  se <- strsplit(header_info, ":")[[1]][2]
  
  return(se)
  
}

logp_retriever <- function(header_info){
  
  logp <- strsplit(header_info, ":")[[1]][3]
  
  return(logp)
  
}

EAF_retriever <- function(header_info){
  
  EAF <- strsplit(header_info, ":")[[1]][4]
  
  return(EAF)
  
}


rsid_retriever <- function(header_info){
  
  rsid <- strsplit(header_info, ":")[[1]][5]
  
  return(rsid)
  
}

##############
#Loading data#
##############

body_fat_data<-read.table("/RAW_DATA/ukb-b-8909.vcf.gz", stringsAsFactors = FALSE)

#We are also gonna load the genome-wide significant SNPs to be able to see if the data
#matches and check whether we have weird results or not with this format.

body_fat_check <- extract_instruments(outcomes = "ukb-b-8909", p1 = 0.00000005, clump = FALSE) #50574

#############
#Small check#
#############

body_fat_check[1,]

#pos.exposure se.exposure samplesize.exposure beta.exposure chr.exposure
#1    210344884  0.00197042              454633     0.0137483            1
#pval.exposure id.exposure        SNP effect_allele.exposure
#1   2.99985e-12  ukb-b-8909 rs78508049                      C
#other_allele.exposure eaf.exposure                             exposure
#1                     T     0.187686 Body fat percentage || id:ukb-b-8909
#mr_keep.exposure pval_origin.exposure data_source.exposure
#1             TRUE             reported                  igd

head(body_fat_data) #the RSID is in V3

body_fat_data[which(body_fat_data$V3 == "rs78508049"),]

#V1        V2         V3 V4 V5 V6   V7          V8             V9
#614290  1 210344884 rs78508049  T  C  . PASS AF=0.187686 ES:SE:LP:AF:ID
#V10
#614290 0.0137483:0.00197042:11.5229:0.187686:rs78508049

#It seems that the effect allele is the second one!!
#I just checked in phenoscanner. And exactly. 

#Let's see if I can transform and de transform the p-value:

abs(log10(body_fat_check$pval.exposure[which(body_fat_check$SNP == "rs78508049")]))  #perfect

#Let's play now how to detransform this:

log10(0.00000005) #-7.3

10^-7.30 #AHA.

#The important thing was the minus.
#That is why...
#We might have an issue here. 

log10(0.05) #it is always gonna be negative.

log10(0.85) #always gonna be negative.

#Hence why people often do the -log10.

#######################
#REFORMATTING THE DATA#
#######################

#Now that we know the basics of the data.
#Let's start by making some functions to make it work:

body_fat_data_copy <- body_fat_data

colnames(body_fat_data_copy) <- c("chr.exposure", "pos.exposure", "SNP", "other_allele.exposure", "effect_allele.exposure", "dot", "PASS", "AF", "Headers", "Headers_info")

#Now we are gonna make functions to retrieve the info on headers info:

body_fat_data_copy$beta.exposure <- as.numeric(as.character(unlist(sapply(body_fat_data_copy$Headers_info, beta_retriever))))
body_fat_data_copy$se.exposure <- as.numeric(as.character(unlist(sapply(body_fat_data_copy$Headers_info, se_retriever))))
body_fat_data_copy$logp <- as.numeric(as.character(unlist(sapply(body_fat_data_copy$Headers_info, logp_retriever))))
body_fat_data_copy$eaf.exposure <- as.numeric(as.character(unlist(sapply(body_fat_data_copy$Headers_info, EAF_retriever))))
body_fat_data_copy$rsid <- as.character(unlist(sapply(body_fat_data_copy$Headers_info, rsid_retriever)))

#To convert to p-value I need to transform them again:

body_fat_data_copy$logp_neg <- body_fat_data_copy$logp*(-1)

#And thus...

body_fat_data_copy$pval.exposure <- 10^(body_fat_data_copy$logp_neg)

############
#FINISHED!!#
############

####################################################################
#Now we start the curation process, let's do some checks beforehand#
####################################################################

#1. Are they OR or logODDs??

summary(body_fat_data_copy$beta.exposure) #No 1s, so no Odd Ratios. Moreover we have SE so...

#2. There is a filter with a dot. What is that?

head(body_fat_data_copy$dot)
tail(body_fat_data_copy$dot)

length(which(body_fat_data_copy$dot == ".")) #all of them.

#something from the vcf file I guess.

#.3 We have a column that says: PASS.

length(which(body_fat_data_copy$PASS == "PASS")) #all of them.

#That means that it was a way to retrieve the good SNPs for the authors.

#4. Let's check again if the info is correct:

body_fat_check[1,]

#pos.exposure se.exposure samplesize.exposure beta.exposure chr.exposure
#1    210344884  0.00197042              454633     0.0137483            1
#pval.exposure id.exposure        SNP effect_allele.exposure
#1   2.99985e-12  ukb-b-8909 rs78508049                      C
#other_allele.exposure eaf.exposure                             exposure
#1                     T     0.187686 Body fat percentage || id:ukb-b-8909
#mr_keep.exposure pval_origin.exposure data_source.exposure
#1             TRUE             reported                  igd

head(body_fat_data) #the RSID is in V3

body_fat_data_copy[which(body_fat_data_copy$SNP == "rs78508049"),]

#chr.exposure pos.exposure        SNP other_allele.exposure
#614290            1    210344884 rs78508049                     T
#effect_allele.exposure dot PASS          AF        Headers
#614290                      C   . PASS AF=0.187686 ES:SE:LP:AF:ID
#Headers_info beta.exposure se.exposure
#614290 0.0137483:0.00197042:11.5229:0.187686:rs78508049     0.0137483  0.00197042
#logp eaf.exposure       rsid logp_neg pval.exposure
#614290 11.5229     0.187686 rs78508049 -11.5229  2.999853e-12

###############
#Perfect match#
###############

#Let's try another one:

body_fat_check$SNP[2]

#"rs12133857"

body_fat_check$beta.exposure[2]

#-0.0131399

body_fat_check$effect_allele.exposure[2]

#"T"

body_fat_check$pval.exposure[2]

#3.50002e-09

#And now we check:

body_fat_data_copy[which(body_fat_data_copy$SNP == "rs12133857"),]

##################
#FREAKING PERFECT#
##################

##########
#CURATION#
##########

#We don't have INFO, so we cannot rely on that. 
#1. But we can remove those variants that are MAF < 0.01

body_fat_data_maf <- body_fat_data_copy[which(body_fat_data_copy$eaf.exposure > 0.01),]
body_fat_data_maf <- body_fat_data_maf[which(body_fat_data_maf$eaf.exposure < 0.99),]

#We go from 9.2M to 7.7M. Good enough!!

#2. Now we remove the MHC region:

body_fat_data_maf_mhc <- body_fat_data_maf[which(as.numeric(body_fat_data_maf$chr.exposure) == 6 & as.numeric(body_fat_data_maf$pos.exposure) >= 26000000 & as.numeric(body_fat_data_maf$pos.exposure) <= 34000000),]

summary(as.numeric(body_fat_data_maf_mhc$chr.exposure)) #perfect.
summary(as.numeric(body_fat_data_maf_mhc$pos.exposure)) #perfect.

#Now let's check if we had any interesting variants there:

summary(as.numeric(body_fat_data_maf_mhc$pval.exposure)) #we actually do. But we said we should remove it so...

body_fat_data_maf_no_mhc <- body_fat_data_maf[which(!(body_fat_data_maf$SNP%in%body_fat_data_maf_mhc$SNP)),]

#3. We should have done that before. But let's check whether we are in build 37. I already know that it is, hence why I just check it here.

head(body_fat_data_maf_no_mhc)

#Perfect.

#4. Let's check the RSIDs:

body_fat_data_maf_no_mhc <- body_fat_data_maf_no_mhc[order(body_fat_data_maf_no_mhc$SNP),]

head(body_fat_data_maf_no_mhc)

#We know that some of them are not in the correct format. 
#So we will try to obtain them.

##########
#CAREFUL!#
##########

#We do have weird SNPs. All the freaking exposure are the same, man.
#Alright, let's check which ones have RSIDs and which ones do not.

rsid_index <- which(str_detect(body_fat_data_maf_no_mhc$SNP, "rs") == TRUE)
no_rsid_index <- which(str_detect(body_fat_data_maf_no_mhc$SNP, "rs") == FALSE) #we have 531

#We have 530.

#In any case, we do not have RSID for these ones.
#But that is not as important.
#We will match with chr_pos. We are gonna check a couple of more no RSID SNPs in dbSNP.
#Just to see what the hell is going on and why there are no RSIDs on them.

body_fat_data_RSID <- body_fat_data_maf_no_mhc[rsid_index,]

#Let's check that this has worked:

head(body_fat_data_RSID) #perfect.

#Now let's do the other one:

body_fat_data_NO_RSID <- body_fat_data_maf_no_mhc[no_rsid_index,]

######################
#Let's check that out#
######################

View(body_fat_data_NO_RSID) #The code detected them properly, as expected.

#Let's get those RSIDs!

chr <- paste("chr", body_fat_data_NO_RSID$chr.exposure, sep = "")
body_fat_data_NO_RSID$chr_pos <- paste(chr, body_fat_data_NO_RSID$pos.exposure, sep = ":")

head(body_fat_data_NO_RSID, 150) #the first one is found in build 37!
tail(body_fat_data_NO_RSID) #the first one is found in build 37!

#Let's see in dbsnp... 
#It seems that we can trust the chr_pos!

body_fat_data_NO_RSID[150,] #ONLY FOUND IN DBSNP FOR BUILD 37. FINAL EVIDENCE MAN.

#We can use the chr_pos to match the data, so no issues here.

######################################################################################
#Thus, let's save the data and match as much as possible with chromosome and position#
######################################################################################

#In the whole process we will be very careful if one of the SNPs RSID is not there.
#Because if that is the case we will have a problem when clumping.

fwrite(body_fat_data_maf_no_mhc, "~/CURATED_DATA/body_fat_UKBB_Combined_Curated_FULL.txt")
