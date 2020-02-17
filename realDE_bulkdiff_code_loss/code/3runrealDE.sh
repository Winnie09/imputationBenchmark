#!/bin/bash -l
#SBATCH --partition=lrgmem
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH -A hji7

ml R
Rscript 03_get_wilcox_fdr.R $1

