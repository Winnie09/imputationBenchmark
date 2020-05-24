bk = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337.count.txt', as.is=T, sep='\t')
tb = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_anno.txt',sep=',')
suppressMessages(library(org.Hs.eg.db))
genename <- select(org.Hs.eg.db, key=as.character(bk$Entrez.Gene.IDs),columns=c("SYMBOL"),keytype="ENTREZID")$SYMBOL
bk = bk[,-1]
bk = as.matrix(bk)
rownames(bk) = genename
bk <- bk[!is.na(row.names(bk)),]
source('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/resource/function.R')
mat = process10x_rmDupGenes(bk)
colnames(mat) = tb[match(sapply(colnames(mat), function(i) paste0(strsplit(i,'_')[[1]][1:2],collapse='_')), tb[,'Sample.name']), 'characteristics.cell.line']
colnames(mat) = paste0(colnames(mat),'_', 1:ncol(mat))
saveRDS(mat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count.rds')

