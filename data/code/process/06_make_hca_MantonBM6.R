library(Matrix)
set.seed(12345)
filelist <- paste0(list.files('/home-4/zji4@jhu.edu/scratch/scdata/HCA/data/immune_census/align/',pattern = '^MantonBM6',full.names = T),'/outs/filtered_gene_bc_matrices/hg19')
smat <- sapply(1:length(filelist),function(id) {
  d <- filelist[id]
  genefile <- paste0(d,"/genes.tsv")
  barcodefile <- paste0(d,"/barcodes.tsv")
  matrixfile <- paste0(d,"/matrix.mtx")
  gn <- read.table(genefile,as.is=T)
  cn <- readLines(barcodefile)
  expr <- as.matrix(readMM(matrixfile))
  row.names(expr) <- paste0(gn[,1],':',gn[,2])
  colnames(expr) <- cn
  expr
})

gene <- unique(unlist(lapply(smat,row.names)))
cell <- unique(unlist(lapply(smat,colnames)))
mat <- matrix(0,nrow=length(gene),ncol=length(cell),dimnames = list(gene,cell))
for (i in smat) {
  mat[row.names(i),colnames(i)] <- i
}

tb = mat
tb <- tb[rowSums(tb) > 0,]
gn = sub('.*:','',rownames(tb))
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
mat = round(tb) ## [1:21180, 1:6941]

## retain cells with at least 500 detected genes
mat = mat[,which(colSums(mat>0) >=500)] ## [1:21180, 1:6939]

## remove low expr genes
mat = mat[which(rowMeans(mat[,] > 0) >=0.01),]  ## [1:12112, 1:6939]
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/MantonBM6/'
saveRDS(mat,paste0(rdir,'genebycell.rds'))
write.table(mat, paste0(rdir,'genebycell.csv'),sep=",",quote=F)
tmat = t(mat)
saveRDS(tmat,paste0(rdir,'cellbygene.rds'))
write.table(tmat,paste0(rdir,'cellbygene.csv'),sep=",",quote=F)

suppressMessages(library(scran))
sce <- SingleCellExperiment(list(counts=mat))
if (ncol(mat) < 21){
      sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5),sizes=c(5,10,15,20))
} else {
      sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5))  
}
sf <- sizeFactors(sce)
normmat <- sweep(mat,2,sf,'/')
normmat <- log2(normmat + 1)
saveRDS(normmat, paste0(rdir,'norm_genebycell.rds'))
write.table(normmat,file=paste0(rdir,'norm_genebycell.csv'),sep=",",quote=F)
tnormmat = t(normmat)
saveRDS(tnormmat, paste0(rdir,'norm_cellbygene.rds'))
write.table(tnormmat,paste0(rdir,'norm_cellbygene.csv'),sep=",",quote=F)
