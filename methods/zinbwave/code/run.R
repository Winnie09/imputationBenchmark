suppressMessages(library(zinbwave))
data <- readRDS(paste0(commandArgs(trailingOnly = T)[1]))
data <- SummarizedExperiment(list(counts=data))
res <- zinbwave(data, K=2, epsilon=1000, normalizedValues=TRUE,residuals = F,BPPARAM=SerialParam())
res <- assays(res)$normalizedValues
res <- log2(exp(res) + 1)
saveRDS(res,paste0(commandArgs(trailingOnly = T)[2]))

