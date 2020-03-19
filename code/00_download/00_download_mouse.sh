#!/bin/bash
#SBATCH --job-name=m_download
#SBATCH --output=m_download.out
#SBATCH --error=m_download.err
#SBATCH --time=6:00:00 
#SBATCH --mem=15000
#SBATCH --partition=rbaltman

# -- FOR rnaseq -- #
#wget https://data-refinery-s3-compendia-circleci-prod.s3.amazonaws.com/MUS_MUSCULUS_1_1574233541.zip

# -- unzip the data -- #
unzip MUS_MUSCULUS_1_1574233541.zip
