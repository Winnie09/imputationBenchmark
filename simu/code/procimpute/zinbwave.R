allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/zinbwave'))
getf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/zinbwave/'))
runf <- setdiff(allf, getf)
res <- lapply(runf, function(f){
      sexpr <-  readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/zinbwave/',f, '.rds'))
      saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/zinbwave/',f,'.rds'))
})

