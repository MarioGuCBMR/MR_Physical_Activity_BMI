# CURATING DATA FOR CAUSE ANALYSIS

The codes found in this folder clean the BFP and WCadjBMI data.

## Data needed:

In /RAW_DATA/ you should have the BFP and WCadjBMI data:

WCadjBMI GWAS summary statistics from Shungin et al 2015 can be found here: https://portals.broadinstitute.org/collaboration/giant/images/7/73/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz
BFP GWAS summary statistics from Elsworth et al 2018 can be found here: https://gwas.mrcieu.ac.uk/datasets/ukb-b-8909/

Also you should have BMI and WHR data from Pulit et al 2019 and Shungin et al 2015, since they are used to clean the data:

BMI GWAS summary statistics from Pulit et al can be found here: https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1
WHR GWAS summary statistics from Shungin et al can be found here: https://portals.broadinstitute.org/collaboration/giant/images/5/54/GIANT_2015_WHR_COMBINED_EUR.txt.gz

## Other tools used:

In parts of the code we have used SNPNexus tool (https://www.snp-nexus.org/v4/)

To obtain the same results follow this protocol:

1) Use build 37/hg19 to query the rsids
2) In population data select 1000G: European.
3) From HapMap select CEU.

Then press enter and download the coordinates in txt format.
