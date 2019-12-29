#!/bin/bash -l
#SBATCH --partition=unlimited
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=2
ml R
Rscript 02_cluster_eachMethod_redo.R $1 $2 

