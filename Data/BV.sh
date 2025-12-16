#!/bin/bash

#SBATCH -J BV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH -o out/BV.%j.out
#SBATCH -e out/BV.%j.err
#SBATCH --time=7-0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jingma@fredhutch.org

ml for
Rscript "BV.R" Gauss01 0 cl 0 0

