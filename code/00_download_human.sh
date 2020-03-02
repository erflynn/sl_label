#!/bin/bash
#SBATCH --job-name=download
#SBATCH --output=logs/download.out
#SBATCH --error=logs/download.err
#SBATCH --time=6:00:00 
#SBATCH --mem=15000
#SBATCH --partition=rbaltman



#wget https://data-refinery-s3-compendia-circleci-prod.s3.amazonaws.com/kfc8oetnjca248dx4lpr805a_ac0d90a1-25c0-4617-838f-212bede63702_compendia.zip

cd humans/
unzip human_compendia.zip
