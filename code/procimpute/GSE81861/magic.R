library(data.table)
sexpr <-  t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/GSE81861/magic/GSE81861_Cell_Line_COUNT.csv'),data.table = F)))
d = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/genebycell.rds')    
sexpr[sexpr<0] <- 0
row.names(sexpr) <- row.names(d)
colnames(sexpr) <- colnames(d)
saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/magic/GSE81861_Cell_Line_COUNT.rds"))
