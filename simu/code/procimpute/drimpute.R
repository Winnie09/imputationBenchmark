library(data.table)
allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/drimpute'))
getf <- sub('.rds','', list.files("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/drimpute/"))
runf <- setdiff(allf, getf)
res <- lapply(runf,function(f) {
  print(f)
  sexpr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/drimpute/',f,'.rds'))
  saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/drimpute/",f,'.rds'))
})

