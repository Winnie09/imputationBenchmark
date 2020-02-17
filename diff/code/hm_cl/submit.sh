#!/bin/bash -l
#SBATCH --partition=lrgmem
#SBATCH --time=36:00:00
#SBATCH --ntasks-per-node=12
ml R
Rscript raw.R
