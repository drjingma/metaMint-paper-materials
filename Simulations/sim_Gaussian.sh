#!/bin/bash

#SBATCH -J sim_Gaussian
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o out/sim_Gaussian.%j.out
#SBATCH -e out/sim_Gaussian.%j.err
#SBATCH --time=2-0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jingma@fredhutch.org
#SBATCH --array=1-100

ml fhR
Rscript "sim_Gaussian.R" ${SLURM_ARRAY_TASK_ID} 100 150 200 3 3 1 1 1 3 0.1 0.8 Gauss Gauss01 3

