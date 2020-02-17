#!/bin/bash -l
#SBATCH --time=10000:0:0
#SBATCH --partition=unlimited
#SBATCH --cpus=24

sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/viper.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/result/impute/viper/$1.rds
