library(Matrix)
process10x_rmDupGenes <- function(genebycellmat){
      tb = genebycellmat
      tb <- tb[rowSums(tb) > 0,]
      gn = rownames(tb)  
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
      tb = round(tb)
}

tb1 = as.matrix(readMM('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/293T/filtered_matrices_mex/hg19/matrix.mtx'))
g1 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/293T/filtered_matrices_mex/hg19/genes.tsv')
barc1 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/293T/filtered_matrices_mex/hg19/barcodes.tsv')
barc1 = paste0('293T_', barc1[,1])
rownames(tb1) = g1[,2]
colnames(tb1) = barc1
libsize1 = colSums(tb1)/1e3
tmp1 = process10x_rmDupGenes(tb1)

tb2 = as.matrix(readMM('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/jurkat/filtered_matrices_mex/hg19/matrix.mtx'))
g2 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/jurkat/filtered_matrices_mex/hg19/genes.tsv')
barc2 = read.table('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/10x/cellline/jurkat/filtered_matrices_mex/hg19/barcodes.tsv')
barc2 = paste0('jurkat_',barc2[,1])
rownames(tb2) = g2[,2]
colnames(tb2) = barc2
libsize2 = colSums(tb2)/1e3
tmp2 = process10x_rmDupGenes(tb2)

## cbind
allg = union(rownames(tmp1),rownames(tmp2))
mat = matrix(0,nrow=length(allg),ncol=sum(ncol(tmp1),ncol(tmp2)))
dimnames(mat) = list(allg, c(colnames(tmp1),colnames(tmp2)))
mat[rownames(tmp1),colnames(tmp1)] = tmp1
mat[rownames(tmp2),colnames(tmp2)] = tmp2

## remove low expr genes
ct = sub('_.*','',colnames(mat))
res <- sapply(unique(ct),function(i){
      kid = which(rowMeans(mat[,ct==i] > 0) >=0.1)
})
names(res) = NULL
kid = unique(names(unlist(res)))
mat = mat[kid,]

mat = round(mat)
saveRDS(mat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/genebycell.rds')
write.table(mat, '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/genebycell.csv',sep=",",quote=F)
tmat = t(mat)
saveRDS(tmat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/cellbygene.rds')
write.table(tmat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/cellbygene.csv',sep=",",quote=F)

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
saveRDS(normmat, '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/norm_genebycell.rds')
write.table(normmat,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/norm_genebycell.csv',sep=",",quote=F)
tnormmat = t(normmat)
saveRDS(tnormmat, '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/norm_cellbygene.rds')
write.table(tnormmat,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/norm_cellbygene.csv',sep=",",quote=F)
