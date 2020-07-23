library(data.table)
allf <- sub('.mat','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/dca'))
getf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/dca'))
runf <- setdiff(allf, getf)
res <- lapply(runf,function(f) {
  if (!f%in%getf){
    print(f)
    sexpr <- fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/dca/',f,'/mean.tsv'),data.table=F)
    row.names(sexpr) <- sexpr[,1]
    sexpr <- as.matrix(sexpr[,-1])
    sexpr = log2(sexpr + 1)
    saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/dca/",f,'.rds'))
  }
})

