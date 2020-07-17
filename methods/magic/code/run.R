# ml python/2.7
library(Rmagic) ## input should be col genes
library(ggplot2)
data = readRDS(commandArgs(trailingOnly = T)[1])
#data <- library.size.normalize(data)
#data <- log2(data+1)
MAGIC_data <- magic(data)
saveRDS(MAGIC_data,commandArgs(trailingOnly = T)[2])
