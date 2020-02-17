#!/bin/bash -l
#SBATCH --time=10:0:0
#SBATCH --mem=80G
#SBATCH --partition=lrgmem
#SBATCH -A hji7
ml R
Rscript 04_get_wilcox_fdr.R $1
