d <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/meta/transfer.csv',as.is=T,header=F)
d <- unique(d[,2])
d <- unlist(sapply(d,function(i) strsplit(i,';')[[1]]))
library(data.table)
tab <- fread('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/meta/metadata.tsv',data.table=F)
tab <- tab[tab[,3]=='gene quantifications' & tab$Assembly=='mm10' & tab[,'File Status']=='released',]
tab <- tab[tab[,4] %in% d,]
url <- tab[,'File download URL']
writeLines(paste0('wget ',url),'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/code/download.sh')
