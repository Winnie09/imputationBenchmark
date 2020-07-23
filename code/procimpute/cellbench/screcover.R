allf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/cellbench/screcover'))
getf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/screcover/'))
runf = setdiff(allf,getf)
runf = c("sc_10x", "sc_celseq2","sc_dropseq","sc_10x_5cl","sc_celseq2_5cl_p1","sc_celseq2_5cl_p2","sc_celseq2_5cl_p3")
res <- lapply(runf, function(f){
  print(f)  
sexpr = read.csv(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/impute/cellbench/screcover/',f, '/scRecover+scImpute.csv'), as.is = T, row.names = 1)
  sexpr = as.matrix(sexpr)
  if (grepl('mix',f)){
    colnames(sexpr) = sapply(colnames(sexpr), function(i) {
      tmp <- strsplit(i,'\\.')[[1]]
      paste0(tmp[1],':',paste(tmp[-c(1,length(tmp))],collapse = '.'),':',tmp[length(tmp)])
    },USE.NAMES = F)  
  } else {
    colnames(sexpr) = sapply(colnames(sexpr), function(i) {
      sub('\\.',':',i)
    },USE.NAMES = F)  
    
  }
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
  saveRDS(normmat,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/screcover/',f,'.rds'))
})
