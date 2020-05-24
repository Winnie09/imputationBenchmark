clumethod = as.character(commandArgs(trailingOnly = T)[[1]])
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/pbmc/',clumethod,'/medianSil/'))
d <- sapply(sub('.rds','',allf), function(f){
  tmp <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/pbmc/',clumethod,'/medianSil/',f,'.rds'))
  if (!is.null(ncol(tmp))){
    tmp
  } else {
    NA
  }
})
d = do.call(rbind,d)
saveRDS(d, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/',clumethod,'_pbmc_allstatistics.rds'))
