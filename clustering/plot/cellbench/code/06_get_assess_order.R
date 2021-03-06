d1 = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/kmeans_cellbench_statistics_scaled_summary.rds'))
# v1 = rowMeans(d1)   ## GB reviewer recommends to average first 3 and then average with the medianSil, 06/22/20
d1_tmp = cbind(mean1 = rowMeans(d1[,1:3]), medianSil = d1[,4])
v1 = rowMeans(d1_tmp)
saveRDS(v1,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/clustering_kmeans_cellbench.rds')
d2 = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/louvein_cellbench_statistics_scaled_summary.rds'))
d2[which(d2<0)]<-0.01 #### 
# v2 = rowMeans(d2)  ## GB reviewer recommends to average first 3 and then average with the medianSil,06/22/20
d2_tmp = cbind(mean1 = rowMeans(d2[,1:3]), medianSil = d2[,4])
v2 = rowMeans(d2_tmp)

saveRDS(v2,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/clustering_louvein_cellbench.rds')

o1 = names(sort(v1))
o2 = names(sort(v2))
saveRDS(o1,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/clustering_kmeans_cellbench.rds')
saveRDS(o2,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/rank/clustering_louvein_cellbench.rds')


