## input: d(matrix), k(number of clusters)
## output: clu(vector)
suppressMessages(library(scran))
suppressMessages(library(igraph))
mtd = as.character(commandArgs(trailingOnly = T)[[1]])
print(mtd)
allf = c("cellmix1","cellmix2","cellmix3","cellmix4","sc_10x","sc_celseq2","sc_dropseq","sc_10x_5cl", "sc_celseq2_5cl_p1", "sc_celseq2_5cl_p2","sc_celseq2_5cl_p3")
df = data.frame(allf = allf, k = c(rep(3,7), rep(5,4)))
getf <- sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd)))
runf <- intersect(allf,getf)
res <- sapply(runf, function(f){
  set.seed(12345)
  print(f)
  d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/',f,'.rds'))
  if (!grepl('latent',mtd)){
    rsd = apply(d,1,sd)
    rm = rowMeans(d)
    cv = rsd/rm
    d = d[cv >= median(cv,na.rm = T), ]  
    d = prcomp(t(d),scale. = T)$x[,1:10] 
  } else {
    d = t(d)  
  }
  ct = sub('.*:','',rownames(d))
  k = as.numeric(df[which(allf==f),'k'])
  graph = buildSNNGraph(d, transposed=T,k=10,d=NA)
  res = cluster_louvain(graph)$membership
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
  dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/cellbench/louvein/',mtd), showWarnings = F)
  saveRDS(clu,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/cellbench/louvein/',mtd,'/',f,'.rds'))
  return(0)
})
