allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/knnsmooth'))
res <- lapply(allf, function(f){
      sexpr <-  readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/impute/knnsmooth/',f,'.rds'))
      suppressMessages(library(scran))
      sce <- SingleCellExperiment(list(counts=sexpr))
      if (ncol(sexpr) < 21){
            sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5),sizes=c(5,10,15,20))
      } else {
            sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5))
      }
      sf <- sizeFactors(sce)
      normmat <- sweep(sexpr,2,sf,'/')
      normmat <- log2(normmat + 1)
      saveRDS(normmat,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/knnsmooth/',f,'.rds'))
})


