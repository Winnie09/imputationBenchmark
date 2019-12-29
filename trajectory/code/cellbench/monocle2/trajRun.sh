#!/bin/bash -l
#SBATCH --partition=lrgmem
#SBATCH --time=4:00:00
#SBATCH --mem=110G
#SBATCH -A hji7

ml R
Rscript  01_get_cor_ov.R $1

