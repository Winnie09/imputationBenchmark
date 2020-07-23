allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/10xcellline/screcover'))
res <- lapply(allf, function(f){
  sexpr = read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/10xcellline/screcover/',f, '/scRecover+scImpute.csv'), as.is = T, row.names = 1)
  sexpr = as.matrix(sexpr)
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
  saveRDS(normmat,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/10xcellline/screcover/',f,'.rds'))
})

