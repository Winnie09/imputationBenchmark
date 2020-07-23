setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
allf <- list.files('./simu/result/impute/saverx')
getf <- list.files("./simu/result/procimpute/saverx/")
runf <- setdiff(allf, getf)
res <- lapply(runf,function(f) {
  print(f)
  resf <- setdiff(list.files(paste0('./simu/result/impute/saverx/',f)), 'data.csv')
  sexpr <- readRDS(paste0('./simu/result/impute/saverx/',f,'/',resf,'/denoised.rds'))$estimate
  d <- readRDS(paste0('./simu/data/processed/diffNumCell/',sub('.rds','',f), '/genebycell.rds'))
  colnames(sexpr) = colnames(d)
  sexpr <- log2(sexpr + 1)
  dir.create('./simu/result/procimpute/saverx/',showWarnings = F, recursive = T)
  saveRDS(sexpr,file=paste0("./simu/result/procimpute/saverx/",f))
})
