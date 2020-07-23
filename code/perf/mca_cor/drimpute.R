library(data.table)
bexpr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/ct.rds')
allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/drimpute'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/mca/drimpute/',f,'.rds'))
  be <- bexpr[,f]
  intgene <- intersect(row.names(sexpr),names(be))
  apply(sexpr[intgene,],2,cor,be[intgene],method='spearman')
})
names(res) <- allf
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/cor/drimpute.rds')

