source('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/methods/knnsmooth/code/code.R')
data <- readRDS(paste0(commandArgs(trailingOnly = T)[1]))
res <- knn_smoothing(data,k=10)
saveRDS(res, paste0(commandArgs(trailingOnly = T)[2]))

