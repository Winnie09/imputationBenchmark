raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/MantonBM6/norm_genebycell.rds')
ct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/ct/MantonBM6.rds')

d <- raw
v = apply(d,1,sd) 
d = d[v >= median(v), ]
d = prcomp(t(d),scale. = T)$x[,1:10] ## cell by PC  

clu <- kmeans(d,ncol(raw)/100)$cluster
cluct <- sapply(1:max(clu),function(sc) {
  tab <- table(ct[clu==sc])
  if (length(tab) > 0) {
    tab <- tab/sum(clu==sc)
    if (max(tab) >= 0.7) {
      names(which.max(tab))
    } else {
      NA
    }
  } else {
     NA
  }
})

res <- data.frame(cell=colnames(raw),cluster=clu,celltype=ct,stringsAsFactors = F)
saveRDS(res, file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cell_celltype_cluster.rds')
df = data.frame(cluster=1:max(clu), celltype=cluct)
saveRDS(df, file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cluster_celltype.rds')

