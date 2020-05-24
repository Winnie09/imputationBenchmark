library(Matrix)
set.seed(12345)
filelist <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/pbmc/')
filelist <- filelist[!grepl('pbmc',filelist)]
smat <- sapply(1:length(filelist),function(id) {
  print(id)
  d <- filelist[id]
  genefile <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/pbmc/',d,"/hg19/genes.tsv")
  barcodefile <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/pbmc/',d,"/hg19/barcodes.tsv")
  matrixfile <- paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/pbmc/',d,"/hg19/matrix.mtx")
  gn <- read.table(genefile,as.is=T)
  cn <- readLines(barcodefile)
  expr <- as.matrix(readMM(matrixfile))
  row.names(expr) <- paste0(gn[,1],':',gn[,2])
  colnames(expr) <- paste0(sub('_filtered_gene_bc_matrices.tar.gz','',d),':',cn)
  expr[rowSums(expr) > 0,]
})

gene <- unique(unlist(lapply(smat,row.names)))
cell <- unique(unlist(lapply(smat,colnames)))
mat <- matrix(0,nrow=length(gene),ncol=length(cell),dimnames = list(gene,cell))
for (i in smat) {
  mat[row.names(i),colnames(i)] <- i
} ## [1:21952, 1:94655] 

tb = mat
tb <- tb[rowSums(tb) > 0,] ## [1:21952, 1:94655]
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
mat = round(tb)  ##  [1:21903, 1:94655]
# > summary(colSums(mat>0))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   142.0   457.0   545.0   626.2   674.0  3765.0 
# > summary(rowSums(mat>0))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#       1      13     254    2706    1738   94632 
## retain cells with at least 500 detected genes
mat = mat[,which(colSums(mat>0) >=500)] ## [1:21903, 1:59620] 

## retain genes expressed in at least 1% cells
mat = mat[which(rowMeans(mat[,] > 0) >=0.01),] # [1:8209, 1:59620]
rdir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/pbmc/sorted/'
dir.create(rdir,showWarnings = F)
saveRDS(mat,paste0(rdir,'genebycell.rds'))
write.table(mat, paste0(rdir,'genebycell.csv'),sep=",",quote=F)
tmat = t(mat)
saveRDS(tmat,paste0(rdir,'cellbygene.rds'))
write.table(tmat,paste0(rdir,'cellbygene.csv'),sep=",",quote=F)
## <-------------- temporarily change row names 
rownames(tmat) = sub('-','_',rownames(tmat)) ##
write.table(tmat,paste0(rdir,'cellbygene.csv'),sep=",",quote=F) ##
## ------------- >

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
