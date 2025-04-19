#!/bin/bash
#SBATCH --array=1-2
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-1:00
#SBATCH -p test
#SBATCH --mem=10G               
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=NONE
module load gcc R #Load R module
R CMD BATCH --quiet --no-restore --no-save "--args $SLURM_ARRAY_TASK_ID" inference_1000g_SBM.R inference_1000g_SBM_${SLURM_ARRAY_TASK_ID}.out  





