clumethod = as.character(commandArgs(trailingOnly = T)[1])
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/cellbench/',clumethod,'/medianSil/'))
res <- lapply(allf,function(f){
  df = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/cellbench/',clumethod,'/medianSil/',f))
  colnames(df)[6] = 'medianSil'
  df
})
d = do.call(rbind,res)
#saveRDS(d, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/assess/','clustering_',clumethod,'_cellbench_allstatistics.rds'))
saveRDS(d, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/',clumethod,'_cellbench_allstatistics.rds'))







