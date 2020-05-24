maxcorcut <- 0.6
difcorcut <- 0
propcut <- 0.7
#maxcorcut <- as.numeric(commandArgs(trailingOnly = T))[1]
#difcorcut <- as.numeric(commandArgs(trailingOnly = T))[2]
#propcut <- as.numeric(commandArgs(trailingOnly = T))[3]

d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/MantonBM6/norm_genebycell.rds') ## scran log2
data <- read.table("/home-4/zji4@jhu.edu/scratch/scdata/HCA/geobulk/GSE74246/data/GSE74246_RNAseq_All_Counts.txt",header=T) 
row.names(data) <- data[,1]
data <- as.matrix(data[,-1])
data <- data[,grep("^X",colnames(data))]
hg19 = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename_genelength.rds')
tgngl <- tapply((hg19$genelength)/1000,hg19$genename,max)
gngl <- as.vector(tgngl)
names(gngl) <- names(tgngl)
data <- data[row.names(data) %in% names(gngl),]
bcount <- data
data <- data/gngl[row.names(data)]

libsize <- colSums(data)
libsize <- libsize/1e6
data <- sweep(data,2,libsize,"/")
data <- data[rowSums(data) > 0,]
data <- log2(data + 1)

ct <- sapply(colnames(data),function(i) strsplit(i,"\\.")[[1]][2],USE.NAMES = F)
intgene <- intersect(row.names(data),row.names(d))
d <- d[intgene,]
data <- data[intgene,]
bcount <- bcount[intgene,]
saveRDS(data, '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/hca/GSE74246_RNAseq_normalized.rds')

bd <- sapply(unique(ct),function(sct) {
  rowMeans(data[,ct==sct])
})

expid <- expand.grid(1:ncol(bd),1:ncol(bd))
expid <- expid[expid[,1]!=expid[,2],]
gs <- lapply(1:nrow(expid),function(i) {
  names(sort(bd[,expid[i,1]]-bd[,expid[i,2]],decreasing = T))[1:100]
})

bdgs <- t(sapply(gs,function(i) {
  colMeans(bd[i,])
}))
dgs <- t(sapply(gs,function(i) {
  colMeans(d[i,])
}))

dgs <- (dgs-rowMeans(dgs))/apply(dgs,1,sd)
bdgs <- (bdgs-rowMeans(bdgs))/apply(bdgs,1,sd)
d <- apply(dgs,2,rank)
bd <- apply(bdgs,2,rank)
library(matrixStats)
scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- apply(data, 1, sd)
  (data - cm) / csd
}

corfunc <- function(m1,m2) {
  scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)
}

cormat <- corfunc(d,bd)
maxcor <- apply(cormat,1,max)
max2cor <- apply(cormat,1,function(i) sort(i,decreasing = T)[2])
ct <- colnames(cormat)[apply(cormat,1,which.max)]
ct[maxcor < maxcorcut] <- NA
ct[maxcor-max2cor < difcorcut] <- NA
names(ct) <- colnames(d)
dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/ct/',showWarnings = F)
saveRDS(ct,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/ct/MantonBM6.rds')  ## cell type

