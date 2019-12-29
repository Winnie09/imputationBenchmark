#!/bin/bash -l
#SBATCH --partition=lrgmem
#SBATCH --time=4:00:00
#SBATCH --mem=30G
#SBATCH -A hji7

ml R
Rscript 07_get_foldchange.R $1

