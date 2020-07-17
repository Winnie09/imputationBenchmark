library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scVI_latent'))
#getf <- sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scVI_latent'))
#runf <- setdiff(allf, getf)
runf <- allf
res <- sapply(runf,function(f) {
    print(f)
    sexpr <- t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scVI_latent/',f,'.csv'),data.table = F))) 
    saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scVI_latent/",f,'.rds'))  
})


