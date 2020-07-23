library(data.table)
allf <- sub('scimpute_count.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/scimpute',pattern = 'scimpute_count.rds'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <-  readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/scimpute/',f,'scimpute_count.rds'))
saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/mca/scimpute/",f,'.rds'))
})

