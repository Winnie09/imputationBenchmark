#!/bin/bash -l
#SBATCH --time=00:30:00
#SBATCH --partition=shared
#SBATCH --ntasks-per-node=1
#SBATCH -A hji7

sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/mcimpute.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/mcimpute/$1.csv
