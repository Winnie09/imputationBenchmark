library(DrImpute)
library(SummarizedExperiment)
data = readRDS(commandArgs(trailingOnly = T)[1])
data.log = log2(data+1)
set.seed(1)
data.imp <- DrImpute(data.log)
saveRDS(data.imp, commandArgs(trailingOnly = T)[2])
# devtools::install_github('gongx030/scDatasets')
# library(scDatasets)
# library(SummarizedExperiment)
# data(usoskin)
# usoskin

