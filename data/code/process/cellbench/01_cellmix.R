suppressMessages(library(scran))
dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/'
savedir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/'
allf = c('cellmix1.count.csv','cellmix2.count.csv','cellmix3.count.csv','cellmix4.count.csv','cellmix5.count.csv')
res <- sapply(allf, function(f){
  print(f)
  savef = sub('.count.csv','',f)
  meta = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.metadata.csv'),as.is=T))
  cnt = as.matrix(read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/raw/cellbench/csv/',savef,'.count.csv'), as.is=T))
  m =  meta[match(colnames(cnt),row.names(meta)),c('H1975', 'H2228', 'HCC827')]
  m =matrix(as.numeric(m),nrow=nrow(m))
  # colnames(cnt) = paste0('H1975_H2228_HCC827:', apply(t(round(apply(m,1,'/',9),2)),1,paste,collapse='_'))
  colnames(cnt) = paste0('H1975_H2228_HCC827:', apply(m,1,paste,collapse='_'))
  colnames(cnt) = paste0(colnames(cnt),':',1:length(colnames(cnt)))
  mtch = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
  mtch = as.matrix(mtch)
  row.names(cnt) = mtch[match(row.names(cnt), sub('\\..*','',mtch[,'geneid'])),'genename']
  mat = cnt[!is.na(row.names(cnt)),]
  source('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/resource/function.R')
  libsize = colSums(mat)/1e3
  mat = process10x_rmDupGenes(mat)
  tmp = mat
  ## remove low expr genes
  kid = which(rowMeans(mat[,] > 0) >=0.1)
  mat = mat[names(kid),]
  mat = round(mat)
  
  if (!file.exists(paste0(savedir,savef))){
    system(paste0('mkdir ',savedir,savef))
  }
  
  saveRDS(mat,paste0(savedir,savef,'/genebycell.rds'))
  write.table(mat, paste0(savedir,savef,'/genebycell.csv'),sep=",",quote=F)
  tmat = t(mat)
  saveRDS(tmat,paste0(savedir,savef,'/cellbygene.rds'))
  write.table(tmat,paste0(savedir,savef,'/cellbygene.csv'),sep=",",quote=F)
  ########### scran
  sce <- SingleCellExperiment(list(counts=mat))
  if (ncol(mat) < 21){
    sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5),sizes=c(5,10,15,20))
  } else {
    sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5))
  }
  sf <- sizeFactors(sce)
  normmat <- sweep(mat,2,sf,'/')
  #############
  # normmat <- sweep(mat,2,libsize,'/') # CPM
  normmat <- log2(normmat + 1)
  saveRDS(normmat, paste0(savedir,savef,'/norm_genebycell.rds'))
  write.table(normmat,file=paste0(savedir,savef,'/norm_genebycell.csv'),sep=",",quote=F)
  tnormmat = t(normmat)
  saveRDS(tnormmat, paste0(savedir,savef,'/norm_cellbygene.rds'))
  write.table(tnormmat,paste0(savedir,savef,'/norm_cellbygene.csv'),sep=",",quote=F)  
  NULL
})
