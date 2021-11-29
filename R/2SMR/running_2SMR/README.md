# Running 2SMR analysis:

Important: to run these analysis you need to have the curated files in the CURATED_DATA folder.
You can obtain those by following the instructions in /R/curating_data_4_2SMR.

**Please add the path where your /CURATED_DATA folder is.**

## Additional files that you might need:

**In the code you will need to add your own path to this data.**

The LD-clumping is performed with ieugwasr::ld_clump_local with the European samples of the 1000G reference panel from gwasglue that can be downloaded here: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
Additionally, you will need to plink binaries. You can check ld_clump_local info for more information: https://rdrr.io/github/MRCIEU/ieugwasr/man/ld_clump_local.html

## Last remarks:

The code is written so that the final files created are saved in the working directory.
Please use setwd() at the beginning of the code to save the data properly.


