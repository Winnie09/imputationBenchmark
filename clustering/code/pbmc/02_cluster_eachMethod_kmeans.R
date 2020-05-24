## input: d(matrix), k(number of clusters)
## output: clu(vector)
suppressMessages(library(scran))
suppressMessages(library(igraph))
mtd = as.character(commandArgs(trailingOnly = T)[[1]])
set.seed(12345)
d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/',mtd,'/sorted.rds'))
if (!grepl('latent', mtd)){
  rsd = apply(d,1,sd)
  rm = rowMeans(d)
  cv = rsd/rm
  d = d[cv >= median(cv,na.rm = T), ]  
  d = prcomp(t(d),scale. = T)$x[,1:10] 
} else{
  d = t(d)
}
tryCatch(clu <- kmeans(d, centers = 10)$cluster,error=function(e){})
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/pbmc/kmeans/',mtd), showWarnings = F,recursive = T)
saveRDS(clu,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/pbmc/kmeans/',mtd,'/sorted.rds'))
