d1 = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/kmeans_pbmc_allstatistics_scaled.rds'))
# v1 = tapply(rowMeans(d1[,3:6]), list(d1[,1]), mean, na.rm=T)
d1 = d1[complete.cases(d1), ]
rownames(d1) = d1[,1]
d1 = cbind(mean1 = rowMeans(d1[,3:5]), medianSil = d1[,6])
v1 = rowMeans(d1)
saveRDS(v1,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/clustering_kmeans_pbmc.rds')

d2 = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/louvein_pbmc_allstatistics_scaled.rds'))
# v2 = tapply(rowMeans(d2[,3:6]), list(d2[,1]), mean, na.rm=T)
d2 = d2[complete.cases(d2), ]
rownames(d2) = d2[,1]
d2 = cbind(mean1 = rowMeans(d2[,3:5]), medianSil = d2[,6])
v2 = rowMeans(d2)

saveRDS(v1,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/clustering_kmeans_pbmc.rds')
saveRDS(v2,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/clustering_louvein_pbmc.rds')
o1 = names(sort(v1))
o2 = names(sort(v2))
saveRDS(o1,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/clustering_kmeans_pbmc.rds')
saveRDS(o2,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/clustering_louvein_pbmc.rds')
