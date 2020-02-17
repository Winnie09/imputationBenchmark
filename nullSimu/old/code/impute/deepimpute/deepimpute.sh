#!/bin/bash -l
#SBATCH --time=6:00:00
#SBATCH --partition=lrgmem
#SBATCH --ntasks-per-node=2
#SBATCH -A hji7

ml R
sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/deepimpute.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/deepimpute/$1.csv
