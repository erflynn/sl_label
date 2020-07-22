#!/bin/bash
#SBATCH --job-name=download_sr
#SBATCH --output=download_sr_%A_%a.out
#SBATCH --error=download_sr_%A_%a.err
#SBATCH --time=10:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman


cd sra/
my_var=`cat ../rnaseq_studies.csv`
for val in $my_var; 
do
    echo $val;
    wget "https://www.ebi.ac.uk/ena/data/view/"${val}"&display=xml"; 

done