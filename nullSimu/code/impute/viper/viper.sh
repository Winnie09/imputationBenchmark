#!/bin/bash -l
#SBATCH --time=4:0:0
#SBATCH --mem=110G
#SBATCH --partition=lrgmem
#SBATCH -A hji7

sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/viper.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1  /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/viper/$1.rds
