# Curating data for 2SMR analysis

The 2SMR analysis use BMI data tha can be donwloaded here: https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz

This data however, lacks chromosome and positions and some variants will not match with the physical activity data because they are merged.
Hence, why this section is so important!

## Data needed to curate the data:

1) Data from Pulit et al 2019 (https://zenodo.org/record/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1)
2) PhenoScanner's dbSNP 147 and 1000G data structured as in https://github.com/MarioGuCBMR/local_SNP_query_in_PhenoScanner.
3) List of merged SNPs from ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz.

## Using PhenoScanner data:

I, unfortunately, cannot distribute dbsnp 147 data used to clean the BMI data used in 2SMR analysis, making curating_data_4_2SMR a bit complicated to reproduce.
For this reason, all the data that this codes produces and that is used in the 2SMR pipeline can be obtained directly from the following zenodo link:


