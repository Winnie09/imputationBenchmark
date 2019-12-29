#!/bin/bash -l
#SBATCH --time=20:0:0
#SBATCH --mem=100G
#SBATCH --partition=lrgmem
#SBATCH -A hji7
ml R
Rscript 03_get_mast_res.R $1
