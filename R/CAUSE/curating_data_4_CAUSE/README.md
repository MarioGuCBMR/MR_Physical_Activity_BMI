# CURATING DATA FOR CAUSE ANALYSIS

The codes found in this folder clean the BFP and WCadjBMI data.

## Data needed:

In /RAW_DATA/ you should have the following BFP, WHR and WCadjBMI data:

WCadjBMI GWAS summary statistics from Shungin et al 2015 can be found here: https://portals.broadinstitute.org/collaboration/giant/images/7/73/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz
BFP GWAS summary statistics from Elsworth et al 2018 can be found here: https://gwas.mrcieu.ac.uk/datasets/ukb-b-8909/
WHR GWAS summary statistics from Shungin et al can be found here: https://portals.broadinstitute.org/collaboration/giant/images/5/54/GIANT_2015_WHR_COMBINED_EUR.txt.gz

Also you should have BMI data from Pulit et al 2019 since it is used to clean the data:

BMI GWAS summary statistics from Pulit et al can be found here: https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1

## Other tools used:

In parts of the code we have used SNPNexus tool (https://www.snp-nexus.org/v4/) and liftover.

To obtain the same results for SNPNexus follow this protocol:

1) Use build 37/hg19 to query the rsids
2) In population data select 1000G: European.
3) From HapMap select CEU.

Then press enter and download the coordinates in txt format.

## Final note:

**You need to generate the curated WHR data and have it in /CURATED_DATA to clean properly the WCadjBMI data.**
