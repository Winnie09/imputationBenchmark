sexpr <- read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/10xcellline/pblr/hg19.csv'),as.is=T, header=F)
sexpr = as.matrix(sexpr)
sexpr = log2(sexpr + 1)
d <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/genebycell.rds')
row.names(sexpr) <- row.names(d)
colnames(sexpr) <- colnames(d)
saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/10xcellline/pblr/hg19.rds'))  
