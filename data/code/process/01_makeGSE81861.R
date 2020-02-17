tb = read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/GSE81861/GSE81861_Cell_Line_COUNT.csv',as.is = T)
rownames(tb) = tb[,1]
tb = as.matrix(tb[,-1])
tb = round(tb)
colnames(tb) <- sub('GM12878_B2','GM12878',colnames(tb))
colnames(tb) <- sub('H1_B2','H1',colnames(tb))
keepcell <- unlist(sapply(c('A549','GM12878','IMR90','H1','K562'),function(i) {
      grep(paste0('__',i,'__'),colnames(tb))
}))
tb <- tb[,keepcell]
colnames(tb) <- sapply(colnames(tb),function(sc) {
      sctmp <- strsplit(sc,'__')[[1]]
      paste0(sctmp[2],'_',sctmp[1])
},USE.NAMES = F)
ct = sub('_.*','',colnames(tb))
libsize = colSums(tb)/1000
tb <- tb[rowSums(tb) > 0,]

gn <- sapply(row.names(tb),function(sc) {
      strsplit(sc,'_')[[1]][2]
},USE.NAMES = F)

rs <- rowSums(tb)
kid <- sapply(unique(gn),function(sid) {
      tmp <- which(gn==sid)
      if (length(tmp)==1) {
            tmp
      } else {
            tmp[which.max(rs[tmp])]
      }
})
tb <- tb[kid,]
row.names(tb) <- gn[kid]
tb <- tb[!grepl('^MT-',row.names(tb)),]


res <- sapply(unique(ct),function(i){
      kid = which(rowMeans(tb[,ct==i] > 0) >=0.1)
})
names(res) = NULL
kid = unique(names(unlist(res)))
tb <- tb[kid,] # 27315   362
saveRDS(tb, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/genebycell.rds'))
write.table(tb,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/genebycell.csv'),sep=",",quote=F)
ttb = t(tb)
saveRDS(ttb, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/cellbygene.rds'))
write.table(ttb,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/cellbygene.csv'),sep=",",quote=F)
########## scran
suppressMessages(library(scran))
sce <- SingleCellExperiment(list(counts=tb))
if (ncol(tb) < 21){
      sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5),sizes=c(5,10,15,20))
} else {
      sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5))  
}
sf <- sizeFactors(sce)
normmat <- sweep(tb,2,sf,'/')
normtb <- log2(normmat + 1)
saveRDS(normtb, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/norm_genebycell.rds'))
write.table(normtb,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/norm_genebycell.csv'),sep=",",quote=F)
tnormtb = t(normtb)
saveRDS(tnormtb, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/norm_cellbygene.rds'))
write.table(tnormtb,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/norm_cellbygene.csv'),sep=",",quote=F)

