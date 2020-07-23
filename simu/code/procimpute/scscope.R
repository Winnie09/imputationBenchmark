
library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scscope'))
#getf = sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scscope/'))
#runf <- setdiff(allf,getf)
runf <- allf
runf <- runf[grepl('_0_0',runf)]
res <- lapply(runf, function(f){
    print(f)
    d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f, '/genebycell.rds'))
    sexpr <-  t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scscope/',f,'.csv'),data.table = F)))
    row.names(sexpr) <- row.names(d)
    colnames(sexpr) <- colnames(d)
    saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scscope/',f,'.rds'))  
})

