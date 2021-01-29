#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --output=logs/counts_mat_%A_%a.out
#SBATCH --error=logs/counts_mat_%A_%a.err
#SBATCH --time=2:00:00 
#SBATCH --mem=10000
#SBATCH --partition=rbaltman

prefix=$1


Rscript code/02_expression/02_rnaseq/02c_get_read_counts.R $prefix $SLURM_ARRAY_TASK_ID
