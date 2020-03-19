#!/bin/bash
#SBATCH --job-name=download
#SBATCH --output=download.out
#SBATCH --error=download.err
#SBATCH --time=6:00:00 
#SBATCH --mem=15000
#SBATCH --partition=rbaltman


# -- FOR microarray -- #
#wget https://data-refinery-s3-compendia-circleci-prod.s3.amazonaws.com/kfc8oetnjca248dx4lpr805a_ac0d90a1-25c0-4617-838f-212bede63702_compendia.zip

# -- FOR rnaseq -- #
#wget https://data-refinery-s3-compendia-circleci-prod.s3.amazonaws.com/HOMO_SAPIENS_1_1574170428.zip

# -- unzip the data -- #
unzip HOMO_SAPIENS_1_1574170428.zip
