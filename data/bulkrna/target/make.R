d <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/meta/transfer.csv',header=F,as.is=T)
d <- d[nchar(d[,2]) > 0,1]
saveRDS(d,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/target/res.rds')
writeLines(d,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/target/res.txt')
