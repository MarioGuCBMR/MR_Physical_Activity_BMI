##############
#INTRODUCTION#
##############

#This is a code to clean the GWAS that we have on BMI from Locke et al from the GIANT consortium.
#We are not going to actually clean it per se..., but we are going to make sure that we save it
#in the most proper way to not bias our results.

#Additionally we will clean the data for the physical activity traits.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
memory.limit(size=80000000)

##################
#Loading function#
##################

parse_rsid <- function(rsid){
  #Function to convert Pulit's data (rsid:A1:A2) to rsid.
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

PS_query <- function(rsid_){
  #Function to obtain chromosome and position from build 37 using dbSNP147 data obtained
  #from PhenoScanner.
  
  #Getting the rsid df ready now
  
  rsid_df_original_copy <- as.data.frame(rsid_)
  rsid_df_original <- as.data.frame(rsid_)
  
  rsid_df_original_copy$chr_pos <- NA
  rsid_df_original$chr_pos <- NA
  
  #Let's get the data from our PS folder.
  
  path_ <- "~/Downloads/dbSNP.txt/dbSNP"
  
  indexing <- seq(0, 298, by = 10)
  #indexing <- indexing[seq(1,3)]
  
  #We will need to parse this by partitioning the jobs:
  
  for(i in indexing){
    
    print(i)
    
    my_files <- c()
    
    if(i == 290){
      
      for(j in seq(i, i+8)){
        
        #print(j)
        
        my_files <- c(my_files, paste(path_, j, sep = "")) 
        
      }
      
    } else {
      
      
      for(j in seq(i, i+9)){
        
        #print(j)
        
        my_files <- c(my_files, paste(path_, j, sep = "")) 
        
      }
      
    }
    
    files_list <- as.list(my_files)
    l <- lapply(files_list, fread, sep=",")
    dt <- rbindlist( l )
    
    #Perfect!!
    #Now we check...
    
    dt_ <- dt[which(dt$rsid%in%rsid_),] 
    
    dt_ <- dt_[which(dt_$a1 != "-"),]
    dt_ <- dt_[which(dt_$a2 != "-"),]
    
    allele_vect <- c("A", "G", "C", "T")
    
    dt_ <- dt_[which(dt_$a1%in%allele_vect),]
    dt_ <- dt_[which(dt_$a2%in%allele_vect),]
    
    rsid_gotcha <- rsid_[which(rsid_%in%dt_$rsid)]
    
    rsid_gotcha <- rsid_gotcha[order(match(rsid_gotcha, dt_$rsid))]
    
    if(is_empty(rsid_gotcha)){
      
      print(i)
      
      next
      
    }
    
    #if(i == 90){
    #  
    #  print("Exception, in PhenoScanner the chr16:500040253 is duplicated. I am gonna keep the one in the online version.")
    #  print("To do so I just need to:")
    #  
    #  dt_ <- dt_[which(duplicated(dt_$hg19_coordinates) == FALSE),]
    #  
    #  rsid_gotcha <- rsid_gotcha[order(match(rsid_gotcha, dt_$hg19_coordinates))]
    #  
    #}
    #
    #if(i == 210){
    #  
    #  print("Exception, in PhenoScanner the chr16:500040253 is duplicated. I am gonna keep the one in the online version.")
    #  print("To do so I just need to:")
    #  
    #  dt_ <- dt_[which(duplicated(dt_$hg19_coordinates) == FALSE),]
    #  
    #  rsid_gotcha <- rsid_gotcha[order(match(rsid_gotcha, dt_$hg19_coordinates))]
    #  
    #}
    
    check <- which(rsid_gotcha != dt_$rsid)
    
    if(is_empty(check) == FALSE){
      
      print("Loop number...")
      print(i)
      print("dind't work")
      
    }
    
    rsid_tmp <- rsid_df_original[which(rsid_df_original$rsid_%in%rsid_gotcha),]
    
    rsid_tmp <- rsid_tmp[order(match(rsid_tmp$rsid_, rsid_gotcha)),]
    
    rsid_tmp$chr_pos <- dt_$hg19_coordinates
    
    rsid_df_old <- rsid_df_original[-which(rsid_df_original$rsid_%in%rsid_gotcha),]
    
    rsid_df_original_ <- rbind(rsid_df_old, rsid_tmp)
    
    rsid_df_original <- rsid_df_original_
    
    my_files <- c()
    
  }
  
  return(rsid_df_original)
  
}

##################
#Loading BMI data#
##################

BMI_pulit <- fread("~/RAW_DATA/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
BMI_locke <- fread("~/RAW_DATA/SNP_gwas_mc_merge_nogc.tbl.uniq.gz")

BMI_pulit$SNP <- as.character(unlist(sapply(BMI_pulit$SNP, parse_rsid)))

##############################
#Getting SNPs are safe to use#
##############################

#This is actually really tricky, but what we are going to do is take all those SNPs that are reported with
#chr:pos 37 in Pulit data.
#Then we are going to explore which of those are in need for some change.
#Let's go.

#########################
#A. Making sure of RSIDs#
#########################

BMI_locke <- BMI_locke[order(BMI_locke$SNP),]

###########################################################
#B. Getting all SNPs that are present in both BMI datasets#
###########################################################

BMI_locke_ok <- BMI_locke[which(BMI_locke$SNP%in%BMI_pulit$SNP),] #almost all

BMI_locke_pot_merged <- BMI_locke[which(!(BMI_locke$SNP%in%BMI_pulit$SNP)),] #3727

#First we are gonna make sure of the chromosome and the position of the potentially merged, because
#they are the complicated ones in this case.

rsid_ <- BMI_locke_pot_merged$SNP

ps_data_locke <- PS_query(rsid_)

#No worries about this, we are going to get the chromosomes and the positions of each SNP
#to do the perfect merging with the accelerator data.

#This worked perfectly, but some of them are weird because they have not been referenced properly with the latest
#IDs.

#Hence, let's do this:

ps_data_locke_correct <- ps_data_locke[which(is.na(ps_data_locke$chr_pos) == FALSE),]
ps_data_locke_incorrect <- ps_data_locke[which(is.na(ps_data_locke$chr_pos) == TRUE),]

#And we are going to check them in phenoscanner online since those report the real rsid, baby.

results_1 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(1,100)])
results_2 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(101,200)])
results_3 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(201,300)])
results_4 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(301,400)])
results_5 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(401,500)])
results_6 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(501,600)])
results_7 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(601,700)])
results_8 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(701,800)])
results_9 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(801,900)])
results_10 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(901,1000)])
results_11 <- phenoscanner::phenoscanner(ps_data_locke_incorrect$rsid_[seq(1000,1009)])

ps_query_df <- rbind(results_1$snps, results_2$snps, results_3$snps,
                     results_4$snps, results_5$snps, results_6$snps,
                     results_7$snps, results_8$snps, results_9$snps,
                     results_10$snps, results_11$snps)

ps_query_df_curated <- ps_query_df %>%
  select(snp, rsid, hg19_coordinates)

colnames(ps_query_df_curated) <- c("oldRSID", "newRSID", "chr_pos")

#Now let's merge it with the rest of the PhenoScanner data:

ps_data_locke_correct$oldRSID <- ps_data_locke_correct$rsid_

ps_data_locke_correct_curated <- ps_data_locke_correct %>%
  select(oldRSID, rsid_, chr_pos)

colnames(ps_data_locke_correct_curated) <- c("oldRSID", "newRSID", "chr_pos")

#Let's bind them together to make a dictionary:

missing_snps_dict <- rbind(ps_data_locke_correct_curated, ps_query_df_curated)

#We lose SNPs..., but hey that is what it is. 
#Those are just gonna be kept as they are.
#Unless... 

#They are merged...
#Let's check that out:

ps_data_locke_incorrect_not_in_PS <- ps_data_locke_incorrect[which(!(ps_data_locke_incorrect$rsid_%in%missing_snps_dict$oldRSID)),] #373 are NOT in PS.
ps_data_locke_incorrect_in_PS <- ps_data_locke_incorrect[which(ps_data_locke_incorrect$rsid_%in%missing_snps_dict$oldRSID),] #636 are in PS

######################################################################################################
#C. Now getting those that need to be potentially changed because they are present in the merged list#
######################################################################################################

merged_rs <- fread("~/RAW_DATA/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

locke_incorrect_not_in_ps_need_changed <- ps_data_locke_incorrect_not_in_PS[which(ps_data_locke_incorrect_not_in_PS$rsid_%in%merged_rs$old_rs),] #27/373 need change.
locke_incorrect_not_in_ps_no_change <- ps_data_locke_incorrect_not_in_PS[which(!(ps_data_locke_incorrect_not_in_PS$rsid_%in%merged_rs$old_rs)),] #346/373 do not need change.

swapping_snps <- merged_rs[which(merged_rs$old_rs%in%locke_incorrect_not_in_ps_need_changed$rsid_),]

swapping_snps <- swapping_snps[order(match(swapping_snps$old_rs, locke_incorrect_not_in_ps_need_changed$rsid_)),]

#Check: 

which(swapping_snps$old_rs != locke_incorrect_not_in_ps_need_changed$rsid_)

locke_incorrect_not_in_ps_need_changed$newRSID <- swapping_snps$new_rs

locke_incorrect_not_in_ps_no_change$newRSID <- locke_incorrect_not_in_ps_no_change$rsid_

locke_merged <- rbind(locke_incorrect_not_in_ps_no_change, locke_incorrect_not_in_ps_need_changed)

colnames(locke_merged)

colnames(locke_merged) <- c("oldRSID", "chr_pos", "newRSID")

curated_snps <- rbind(locke_merged, missing_snps_dict)

#PERECT. #Let's check the duplicate:

curated_snps <- curated_snps[order(curated_snps$oldRSID),]

which(duplicated(curated_snps$oldRSID) == TRUE) #we can remove this bad boy.

curated_snps <- curated_snps[which(duplicated(curated_snps$oldRSID) == FALSE),]

#OKAY we have the dictionary ready.

#We have combined:

#1) those that are in MY version of PS locally.
#2) combined those that are not there, but are in the internet version.
#3) we have changed the 27/372 that were not found in either of them, but that have been reported
#as being merged.

fwrite(curated_snps, "~/CURATED_DATA/BMI_locke_dictionary_PS.txt")

#Nowe we are gonna develop the dictionary for those that in Pulit:

Pulit_dictionary <- BMI_pulit[which(BMI_pulit$SNP%in%BMI_locke_ok$SNP),]

Pulit_dictionary$chr_ <- paste("chr", Pulit_dictionary$CHR, sep = "")
Pulit_dictionary$chr_pos <- paste(Pulit_dictionary$chr_, Pulit_dictionary$POS, sep = ":")

Pulit_dictionary_ <- Pulit_dictionary %>%
  select(SNP, chr_pos)

colnames(Pulit_dictionary_) <- c("oldRSID", "chr_pos")

Pulit_dictionary_$newRSID <- Pulit_dictionary_$oldRSID

final_dictionary <- rbind(Pulit_dictionary_, curated_snps)

fwrite(final_dictionary, "~/CURATED_DATA/BMI_locke_FULL_dictionary.txt")

dict_test <- final_dictionary[which(final_dictionary$oldRSID%in%locke_incorrect_not_in_ps_need_changed$rsid_),]
dict_test <- final_dictionary[which(final_dictionary$oldRSID%in%missing_snps_dict$oldRSID),]

######################################
#NICE: recapitulation of what we have#
######################################

#Finally let's save this data:

BMI_locke_dupl <- BMI_locke[order(BMI_locke$p),]
BMI_locke_dupl <- BMI_locke_dupl[which(duplicated(BMI_locke_dupl$SNP) == FALSE),]

dict_match_locke <- final_dictionary[which(final_dictionary$oldRSID%in%BMI_locke_dupl$SNP),]
dict_match_locke <- dict_match_locke[which(duplicated(dict_match_locke$oldRSID) == FALSE),]
dict_match_locke <- dict_match_locke[order(match(dict_match_locke$oldRSID, BMI_locke_dupl$SNP)),]

#Check: 

which(dict_match_locke$oldRSID != BMI_locke_dupl$SNP)

BMI_locke_dupl$chr_pos <- dict_match_locke$chr_pos
BMI_locke_dupl$newRSID <- dict_match_locke$newRSID

fwrite(BMI_locke_dupl, "~/CURATED_DATA/BMI_Locke_Chr_Pos.txt")

###################################################################################
#Now we are gonna check those that are OK because they are present in each PA file#
###################################################################################

#We are just going to do exactly the same as the last step:

#Order by pval.
#Removing duplicates.
#And, finally, saving data with chromosome and position.

#Tomorrow I wake up early and I run this before the meeting with Hermina and Tuomas 
#(which might or might not happen.)

#######################
#Thus, loading PA data#
#######################

sedentary <- fread("~/RAW_DATA/Doherty-2018-NatureComms-sedentary.csv.gz")
moderate <- fread("~/RAW_DATA/Doherty-2018-NatureComms-moderate.csv.gz")
acc425 <- fread("~/RAW_DATA/Acc425_Model1_BOLTLMM_500K.txt.gz")

sedentary_ <- sedentary[order(sedentary$P_BOLT_LMM_INF),]
moderate_ <- moderate[order(moderate$P_BOLT_LMM_INF),]
acc425_ <- acc425[order(acc425$P_BOLT_LMM_INF),]

#Now we remove duplicates by chromosome and position:

sedentary_$chr_ <- paste("chr", sedentary_$CHR, sep = "")
sedentary_$chr_pos <- paste(sedentary_$chr_, sedentary_$BP, sep = ":")

moderate_$chr_ <- paste("chr", moderate_$CHR, sep = "")
moderate_$chr_pos <- paste(moderate_$chr_, moderate_$BP, sep = ":")

acc425_$chr_ <- paste("chr", acc425_$CHR, sep = "")
acc425_$chr_pos <- paste(acc425_$chr_, acc425_$BP, sep = ":")

#Finally we remove duplicates as 2SMR does, but by chromosome and position:

sedentary_dupl <- sedentary_[which(duplicated(sedentary_$chr_pos) == FALSE),]
moderate_dupl <- moderate_[which(duplicated(moderate_$chr_pos) == FALSE),]
acc425_dupl <- acc425_[which(duplicated(acc425_$chr_pos) == FALSE),]

fwrite(sedentary_dupl, "~/CURATED_DATA/Curated_GWAS/Sedentary_doherty.txt")
fwrite(moderate_dupl, "~/CURATED_DATA/Moderate_doherty.txt")
fwrite(acc425_dupl, "~/CURATED_DATA/acc425_klimentidis.txt")



