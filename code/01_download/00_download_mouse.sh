#!/bin/bash
#SBATCH --job-name=m_download
#SBATCH --output=m_download.out
#SBATCH --error=m_download.err
#SBATCH --time=48:00:00 
#SBATCH --mem=15000
#SBATCH --partition=rbaltman

# -- FOR rnaseq -- #
my_f="MUS_MUSCULUS_1_1574233541.zip"
wget https://data-refinery-s3-compendia-circleci-prod.s3.amazonaws.com/${my_f}

# -- unzip the data -- #
unzip ${my_f}
