#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --partition=shared
#SBATCH --mem=110G
#SBATCH -A hji7

sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/pblr.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/pblr/$1.csv
