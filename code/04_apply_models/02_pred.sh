#!/bin/bash
#SBATCH --job-name=pred_%A_%a
#SBATCH --output=logs/pred_%A_%a.out
#SBATCH --error=logs/pred_%A_%a.err
#SBATCH --time=1:00:00 
#SBATCH --mem=10000
#SBATCH --partition=rbaltman

prefix=$1
i=$SLURM_ARRAY_TASK_ID

ml R/3.6.1

Rscript code/04a_pred_sl.R $prefix $i
