library(data.table)
allf <- sub('scimpute_count.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/10xcellline/scimpute',pattern = 'scimpute_count.rds'))
res <- lapply(allf,function(f) {
  print(f)
  sexpr <-  readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/10xcellline/scimpute/',f,'scimpute_count.rds'))
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
  saveRDS(normmat,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/10xcellline/scimpute/",f))
})

