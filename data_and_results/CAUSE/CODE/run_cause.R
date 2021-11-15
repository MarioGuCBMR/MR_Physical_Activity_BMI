##############
#INTRODUCTION#
##############

#This code allows you to:

#1) Simple check on the results and the dataframes used to obtain them.
#2.1) re-run the last step of the CAUSE pipeline with the data available of each combination of traits 
#(it uses the merged_gwas, params and pruned_snps to obtain the results)
#2.2) the same as point 2, but using only the final results and the parameters 
#(no need to use the zenodo files or the pruned_snps.)

#IMPORTANT TO HAVING INSTALLED CAUSE IN THE FOLLOWING WAY:

#devtools::install_version("mixsqp", version = "0.1-97", repos = "http://cran.us.r-project.org")
#devtools::install_version("ashr", version = "2.2-32", repos = "http://cran.us.r-project.org")
#devtools::install_github("jean997/cause")

###################
#Loading libraries#
###################

library(data.table)
library(cause)
library(tidyverse)

##############
#Loading data#
##############

#The merged_df can be obtained in the zenodo link available in the repository README.
#The rest of the files can be obtained in the https://github.com/MarioGuCBMR/MR_Physical_Activity_BMI/tree/main/data_and_results/CAUSE 
#folder. 

#Here is an example:

#merged_df <- readRDS("~/BMI_analysis//Vigorous_BMI/X_Vigorous_BMI.rds") #only available through zenodo because it is too big.
#pruned_snps_df <- read.table("~/BMI_analysis/Vigorous_BMI/Vigorous_BMI_pruned.txt", header = TRUE)
#params <- readRDS("~/BMI_analysis/Vigorous_BMI/params_Vigorous_BMI.rds")
#results <- readRDS("~/BMI_analysis/Vigorous_BMI/RES_Vigorous_BMI.rds")

#NOTE: the pruned snps is a dataframe with the summary statistics of the exposure trait. In this case it contain the summary statistics of Vigorous Physical Activity (acc425) for 1399 pruned variants.
#for this reason:

#pruned_snps <- pruned_snps_df$SNP

#Here you can load yours:

#merged_df <- readRDS("path") #only available through zenodo because it is too big.
pruned_snps_df <- read.table("path", header = TRUE)
params <- readRDS("path")
results <- readRDS("path")

###########################
#STEP 1: check the results#
###########################

#First we check the ELPD results:

results$elpd

#Then, the CAUSE results. The p-values reflects whether the causal model is a better fit.
#Note that the p-values might slightly differ for every run, but they should all be in a similar order.
#However, the median joint posterior probabilities of the casual effect, shared effect and the proportion of shared variants should be the same.

summary(results)

#And with this you can check the sensitivity plots:

plot(res, type="data")

############################
#STEP 2: re-run the results#
############################

#IMPORTANT: this version of this process can only be run with the merged dataframes directly (all variants available for both traits)
#In any case, this first version just shows how it would be run normally.
#If you have not downloaded the zenodo files, go to the next version (faster).

res <- cause(X=merged_df, variants = pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

#Now you can check:

#First we check the ELPD results:

res$elpd

#Then, the CAUSE results. The p-values reflects whether the causal model is a better fit.
#Note that the p-values might slightly differ for every run, but they should all be in a similar order.
#However, the median joint posterior probabilities of the casual effect, shared effect and the proportion of shared variants should be the same.

summary(res)

#And with this you can check the sensitivity plots:

plot(res, type="data")

################################
#STEP 2.2: without zenodo files#
################################

res <- cause(X=results$data, variants = results$data$snp, param_ests = params, qalpha = 1, qbeta = 2)

#Now you can check:

#First we check the ELPD results:

res$elpd

#Then, the CAUSE results. The p-values reflects whether the causal model is a better fit.
#Note that the p-values might slightly differ for every run, but they should all be in a similar order.
#However, the median joint posterior probabilities of the casual effect, shared effect and the proportion of shared variants should be the same.

summary(res)

#And with this you can check the sensitivity plots:

plot(res, type="data")