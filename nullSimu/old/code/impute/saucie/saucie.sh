#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --mem=50G
#SBATCH --partition=shared
#SBATCH -A hji7
mkdir /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/saucie_latent
sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/shell/saucie.sh /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/$1 /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/saucie/$1.csv /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/saucie_latent/$1.csv
