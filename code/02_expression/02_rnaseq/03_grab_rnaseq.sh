#!/bin/bash
#SBATCH --job-name=rnaseq_mat
#SBATCH --output=logs/rnaseq_mat_%A_%a.out
#SBATCH --error=logs/rnaseq_mat_%A_%a.err
#SBATCH --time=2:00:00 
#SBATCH --mem=20000
#SBATCH --partition=rbaltman

prefix=$1
sample_f=$2
outfile=$3

Rscript code/02_expression/02_rnaseq/03_grab_rnaseq_df.R $prefix $sample_f $outfile
