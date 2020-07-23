  sexpr <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/GSE81861/pblr/GSE81861_Cell_Line_COUNT.csv'),as.is=T,header=F)
  sexpr = as.matrix(sexpr)
  sexpr = log2(sexpr + 1)
  d = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/genebycell.rds')
  row.names(sexpr) <- row.names(d)
  colnames(sexpr) <- colnames(d) 
  saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/pblr/GSE81861_Cell_Line_COUNT.rds'))   

