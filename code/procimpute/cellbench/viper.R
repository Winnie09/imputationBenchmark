allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/cellbench/viper'))
res <- lapply(allf, function(f){
      sexpr <-  readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/cellbench/viper/',f, '.rds'))
      saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/viper/',f,'.rds'))
})

