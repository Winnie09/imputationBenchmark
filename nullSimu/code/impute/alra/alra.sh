#!/bin/bash -l
#SBATCH --time=10:0:0
#SBATCH --mem=50G
#SBATCH --partition=skylake
#SBATCH -A hji7

ml R
sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/alra.sh  /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1  /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/alra/$1.rds

