#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --mem=110G
#SBATCH --partition=shared
#SBATCH -A hji7

ml R
sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/drimpute.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/drimpute/$1.rds


