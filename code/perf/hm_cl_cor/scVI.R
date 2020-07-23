library(data.table)
bexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
sexpr <- t(as.matrix(fread('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/GSE81861/scVI/GSE81861_Cell_Line_COUNT.csv',data.table=F)))
d = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/genebycell.rds')
row.names(sexpr) <- row.names(d)

res <- lapply(colnames(bexpr), function(cl){
  be <- bexpr[,cl]
  intgene <- intersect(row.names(sexpr),names(be))
  apply(sexpr[intgene,],2,cor,be[intgene],method='spearman')
})
names(res) <- colnames(bexpr)
saveRDS(res,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/hm_cellline_cor/scVI.rds')

