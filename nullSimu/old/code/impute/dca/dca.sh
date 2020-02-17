#!/bin/bash -l
#SBATCH --time=1:0:0
#SBATCH --partition=shared
#SBATCH --ntasks-per-node=2
#SBATCH -A hji7

sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/dca.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/dca/$1
