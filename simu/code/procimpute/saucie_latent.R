library(data.table)
allf <- sub('.csv','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/saucie_latent'))
#getf <- sub('.rds','', list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/saucie_latent/'))
#runf <- setdiff(allf, getf)
runf <- allf
#sink('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/code/procimpute/saucie_latent_redo.txt')
res <- sapply(runf, function(f){
      sexpr <-  t(as.matrix(fread(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/impute/saucie_latent/',f,'.csv'),data.table = F)))
      saveRDS(sexpr,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/saucie_latent/",f,'.rds'))  
})
#sink()

