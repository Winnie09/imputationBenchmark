#!/bin/bash -l
#SBATCH --time=1:0:0
#SBATCH --mem=50G
#SBATCH --partition=shared

sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/baynorm.sh  /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/baynorm/$1.rds

