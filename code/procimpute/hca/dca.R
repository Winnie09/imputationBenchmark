library(data.table)
allf <- sub('.mat','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/hca/dca'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/hca/dca/',f,'/mean.tsv'),data.table=F)
  row.names(sexpr) <- sexpr[,1]
  sexpr <- as.matrix(sexpr[,-1])
  sexpr = log2(sexpr + 1)
saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/dca/",f,'.rds'))
})

