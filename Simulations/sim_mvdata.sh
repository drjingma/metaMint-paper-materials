#!/bin/bash

#SBATCH -J sim_mvdata
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o out/mvdata.%j.out
#SBATCH -e out/mvdata.%j.err
#SBATCH --time=2-0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jingma@fredhutch.org
#SBATCH --array=1-100

ml fhR
Rscript "sim_mvdata.R" ${SLURM_ARRAY_TASK_ID} 100 49 50 3 3 200 4 1 0 0.1 0.8 ln Gauss01 0

