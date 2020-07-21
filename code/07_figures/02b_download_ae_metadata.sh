#!/bin/bash
#SBATCH --job-name=download_ae
#SBATCH --output=download_ae_%A_%a.out
#SBATCH --error=download_ae_%A_%a.err
#SBATCH --time=2:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman


cd array_express/
my_var=`cat ../non_gse.csv`
for val in $my_var; 
do
    echo $val;
    wget https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/${val};
done