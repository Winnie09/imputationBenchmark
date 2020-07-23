bexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
sexpr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/GSE81861/saver/GSE81861_Cell_Line_COUNT.rds')$estimate
res <- lapply(colnames(bexpr), function(cl){
  be <- bexpr[,cl]
  intgene <- intersect(row.names(sexpr),names(be))
  apply(sexpr[intgene,],2,cor,be[intgene],method='spearman')
})

names(res) <- colnames(bexpr)
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/hm_cellline_cor/saver.rds')

