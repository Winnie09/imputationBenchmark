#!/bin/bash -l
#SBATCH --partition=lrgmem
#SBATCH --time=20:00:00
#SBATCH --mem=100G
ml R
Rscript 02_cluster_eachMethod_louvein.R $1
#Rscript 05_cluster_eachMethod_kmeans.R $1 

