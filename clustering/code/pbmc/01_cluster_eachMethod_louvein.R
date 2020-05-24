## input: d(matrix), k(number of clusters)
## output: clu(vector)
suppressMessages(library(scran))
suppressMessages(library(igraph))
mtd = as.character(commandArgs(trailingOnly = T)[[1]])
print(mtd)
set.seed(12345)
d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/',mtd,'/sorted.rds'))
if (!grepl('latent',mtd)){
  rsd = apply(d,1,sd)
  rm = rowMeans(d)
  cv = rsd/rm
  d = d[cv >= median(cv,na.rm = T), ]  
  d = prcomp(t(d),scale. = T)$x[,1:10] 
} else {
  d = t(d)  
}
ct = sub(':.*','',rownames(d))
graph = buildSNNGraph(d, transposed=T,k=10,d=NA)
res = cluster_louvain(graph)$membership
k = 10
if (max(res) <= 3){
  hclu <- hclust(dist(d))
  clu <- cutree(hclu,k)
} else {
  cc <- aggregate(d,list(res),mean)
  cc <- as.matrix(cc[,-1])
  hclu <- hclust(dist(cc))
  clu <- cutree(hclu,k)
  clu <- clu[res]      
}
names(clu) = row.names(d)
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/pbmc/louvein/',mtd), showWarnings = F)
saveRDS(clu,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/pbmc/louvein/',mtd,'/sorted.rds'))

