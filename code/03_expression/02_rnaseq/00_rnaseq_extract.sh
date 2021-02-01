#!/bin/bash
#SBATCH --job-name=rnaseq_mat
#SBATCH --output=logs/rnaseq_mat_%A_%a.out
#SBATCH --error=logs/rnaseq_mat_%A_%a.err
#SBATCH --time=6:00:00 
#SBATCH --mem=10000
#SBATCH --partition=rbaltman

prefix=$1
#i=$SLURM_ARRAY_TASK_ID


#Rscript code/02e_grab_rnaseq.R $prefix $i

#Rscript code/02g_check_output.R $prefix
Rscript code/03_expression/02_rnaseq/01_quant_to_csv.R $prefix $SLURM_ARRAY_TASK_ID
