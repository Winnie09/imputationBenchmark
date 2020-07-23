library(data.table)
sexpr <-  t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/10xcellline/magic/hg19.csv'),data.table = F)))
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/genebycell.rds')
sexpr[sexpr<0] <- 0
row.names(sexpr) <- row.names(d)
colnames(sexpr) <- colnames(d)
saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/10xcellline/magic/hg19.rds"))

