library(data.table)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
bulk = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count.rds')
g = fread('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/gencode.v19.annotation.gtf',data.table = F)
g <- g[g[,3]=='gene',]
gn <- gsub('\"','',sub(' gene_name ','',sapply(g[,9],function(i) strsplit(i,';')[[1]][5])))
gl <- g[,5]-g[,4]+1
names(gl) <- gn
gl <- gl/1000
bulk <- bulk[row.names(bulk) %in% names(gl),]
bulk <- bulk/gl[row.names(bulk)]
lib <- colSums(bulk)/1e6
bulk <- t(t(bulk)/lib)
bulk <- log2(bulk + 1) ## TPM
colnames(bulk) = sub('_.*','',colnames(bulk))
bulk <- sapply(unique(colnames(bulk)),function(i) rowMeans(bulk[,colnames(bulk)==i]))
saveRDS(bulk,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count_average_replicates.rds')
