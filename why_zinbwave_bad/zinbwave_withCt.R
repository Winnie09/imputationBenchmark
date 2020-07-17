suppressMessages(library(zinbwave))
data <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/genebycell.rds')
ct <- sub('_.*','',colnames(data))
id <- c(which(ct==unique(ct)[1])[1:50],which(ct==unique(ct)[2])[1:50])
data <- data[,id]
ct=ct[id]
data = data[colMeans(data>0)>=0.1, ]
data <- SummarizedExperiment(list(counts=data))
x <- cbind(1,as.numeric(as.factor(ct))-1)
res <- zinbwave(data, X=x, K=2, epsilon=1000, normalizedValues=TRUE,residuals = F,BPPARAM=SerialParam())
#res <- assays(res)$normalizedValues
#res <- log2(exp(res) + 1)
saveRDS(res,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/zinbwave_10x_withCt.rds')




