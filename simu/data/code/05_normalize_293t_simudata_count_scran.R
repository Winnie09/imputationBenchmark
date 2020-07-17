# the data is 'mat'
## /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/code/05_normalize_293t_simudata_count_scran.R
numcell1 = as.numeric(commandArgs(trailingOnly = T)[1])
numcell2 = as.numeric(commandArgs(trailingOnly = T)[2])
probGene = as.numeric(commandArgs(trailingOnly = T)[3])
probRead= as.numeric(commandArgs(trailingOnly = T)[4])
allf = NULL
if (numcell1 <= numcell2){
  if (probRead == 0 | probGene == 0){
    probGene <- probRead <- 0
  } 
  simumat = readRDS(file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/genebycell.rds'))
  suppressMessages(library(scran))
  sce <- SingleCellExperiment(list(counts=simumat))
  if (ncol(simumat) < 21){
    sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5),sizes=c(5,10,15,20))
  } else {
    sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5))
  }
  sf <- sizeFactors(sce)
  normmat <- sweep(simumat,2,sf,'/')
  normmat <- log2(normmat + 1)
  saveRDS(normmat, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/norm_genebycell.rds'))
  write.table(normmat,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/norm_genebycell.csv'),sep=",",quote=F)
  tnormmat = t(normmat)
  saveRDS(tnormmat, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/norm_cellbygene.rds'))
  write.table(tnormmat,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',numcell1,'_',numcell2,'_',probGene,"_",probRead,'/norm_cellbygene.csv'), sep=",",quote=F)
}

