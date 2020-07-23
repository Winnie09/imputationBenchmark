library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scVI'))
getf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scVI'))
runf <- setdiff(allf, getf)
res <- lapply(runf,function(f) {
    print(f)
    d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f, '/genebycell.rds'))
    sexpr <- t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scVI/',f,'.csv'),data.table = F))) 
    sexpr <- log2(sexpr + 1)
    row.names(sexpr) <- row.names(d)
    colnames(sexpr) <- colnames(d)
    saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scVI/",f,'.rds'))  
})


