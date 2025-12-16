#!/bin/bash

#SBATCH -J sim_Gamma
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -o out/sim_Gamma.%j.out
#SBATCH -e out/sim_Gamma.%j.err
#SBATCH --time=2-0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jingma@fredhutch.org
#SBATCH --array=1-100

ml fhR
Rscript "sim_Gamma.R" ${SLURM_ARRAY_TASK_ID} 100 150 200 3 3 1 0.5 4 8 0.1 0.2 Gamma Exp1Gamma 1

