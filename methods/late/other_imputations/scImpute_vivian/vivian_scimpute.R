# install.packages("devtools")
# library(devtools)
# install_github("Vivianstats/scImpute")

# PBMC
library(scImpute)
getwd()

scimpute(count_path='1M_neurons_matrix_subsampled_20k.csv.gz',
         infile = "csv",           # format of input file
         outfile = "csv",          # format of output file
         out_dir = "./",           # full path to output directory
         labeled = FALSE,          # cell type labels not available
         drop_thre = 0.5,          # threshold set on dropout probability
         Kcluster = 2,             # 2 cell subpopulations
         ncores = 10)
