suppressMessages(library(scran))
dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/'
allfd = list.files(dir)
for (fd in allfd){
      simumat <- readRDS(paste0(dir,fd,'/genebycell.rds'))
      sce <- SingleCellExperiment(list(counts=simumat))
      if (ncol(simumat) < 21){
            sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5),sizes=c(5,10,15,20))
      } else {
            sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5))
      }
      sf <- sizeFactors(sce)
      normmat <- sweep(simumat,2,sf,'/')
      normmat <- log2(normmat + 1)
      saveRDS(normmat, paste0(dir,fd,'/norm_genebycell.rds'))
      write.table(normmat,file=paste0(dir,fd,'/norm_genebycell.csv'),sep=",",quote=F)
      tnormmat = t(normmat)
      saveRDS(tnormmat, paste0(dir,fd,'/norm_cellbygene.rds'))
      write.table(tnormmat,paste0(dir,fd,'/norm_cellbygene.csv'), sep=",",quote=F)      
}

