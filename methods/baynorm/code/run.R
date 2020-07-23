library(bayNorm)
data <- readRDS(paste0(commandArgs(trailingOnly = T)[1]))
data <- SummarizedExperiment::SummarizedExperiment(assays=list(Counts=data))
res <- bayNorm(
    Data=data,
    BETA_vec = NULL,
    mode_version=TRUE,
    mean_version = FALSE,S=20
    ,verbose =FALSE,
    parallel = TRUE)
saveRDS(res, paste0(commandArgs(trailingOnly = T)[2]))
