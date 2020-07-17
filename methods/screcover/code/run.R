library(scRecover)
library(BiocParallel)
suppressMessages(library(SingleCellExperiment))
data <- readRDS(paste0(commandArgs(trailingOnly = T)[1]))
scRecover(counts = data, Kcluster = 2, outputDir = paste0(commandArgs(trailingOnly = T)[2],'/'), verbose = FALSE)

