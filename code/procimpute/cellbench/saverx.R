setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result')
allf <- list.files('./impute/cellbench/saverx')
res <- lapply(allf,function(f) {
  print(f)
  resf = setdiff(list.files(paste0('./impute/cellbench/saverx/',f)),'data.csv')
  sexpr <- readRDS(paste0('./impute/cellbench/saverx/',f,'/',resf,'/denoised.rds'))$estimate
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
  sexpr <- log2(sexpr + 1)
  dir.create("./procimpute/cellbench/saverx/", showWarnings = F, recursive = T)
  saveRDS(sexpr,file=paste0("./procimpute/cellbench/saverx/",f))
})

