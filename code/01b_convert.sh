#!/bin/bash
#SBATCH --job-name=gconvert_%A_%a
#SBATCH --output=logs/gconvert_%A_%a.out
#SBATCH --error=logs/gconvert_%A_%a.err
#SBATCH --time=6:00:00 
#SBATCH --mem=40000
#SBATCH --partition=rbaltman

prefix=$1
i=$SLURM_ARRAY_TASK_ID


/scratch/users/erflynn/applications/miniconda3/etc/profile.d/conda.sh;
source activate cmappy;


python3.7 ~/gct_to_gctx.py -f ${prefix}_${i}.gct;
