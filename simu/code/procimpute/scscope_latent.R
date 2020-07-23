
library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scscope_latent'))
#getf = sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scscope_latent/'))
#runf <- setdiff(allf,getf)
runf <- allf
runf <- runf[grepl('_0_0',runf)]
res <- lapply(runf, function(f){
    print(f)
    sexpr <-  t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/scscope_latent/',f,'.csv'),data.table = F)))
    saveRDS(sexpr,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/scscope_latent/',f,'.rds'))  
})

