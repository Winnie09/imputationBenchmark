# [whou10@jhu.edu@bc-login03 procimpute]$ cat scimpute.R
suppressMessages(library(data.table))
allf <- sub('scimpute_count.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scimpute',pattern = 'scimpute_count.rds'))
getf <- list.files("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scimpute/")
runf <- setdiff(allf, getf)
res <- lapply(runf,function(f) {
      print(f)
      sexpr <-  readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scimpute/',f,'scimpute_count.rds'))
      d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',sub('.rds','',f), '/genebycell.rds'))
      suppressMessages(library(scran))
      sce <- SingleCellExperiment(list(counts=d))
      flag <- 1
      if (ncol(d) < 21){
            tryCatch({sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5),sizes=c(5,10,15,20), positive=TRUE);flag <- 0},warning=function(w){},error=function(e){})
      } else {
            tryCatch({sce <- computeSumFactors(sce,BPPARAM=MulticoreParam(workers = 5), positive=TRUE);flag <- 0},warning=function(w){},error=function(e){})
      }
      if (flag==0) {
        sf <- sizeFactors(sce)  
      } else {
        sf <- colSums(d)
        sf <- sf/median(sf)
      }
      normmat <- sweep(sexpr,2,sf,'/')
      normmat <- log2(normmat + 1)
      saveRDS(normmat,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scimpute/",f))
})

