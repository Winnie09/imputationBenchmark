realfdr <- sapply(list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_mast/diff/'),function(method) {
      af <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_mast/diff/',method))
      sapply(af,function(f) {
            d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_mast/diff/',method,'/',f))
            c(approx(x=d[,1],y=d[,2],xout = 0.1)$y)
      })      
})


reportfdr <- sapply(list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_mast/diff/'),function(method) {
      af <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_mast/diff/',method))
      sapply(af,function(f) {
            d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_mast/diff/',method,'/',f))
            c(approx(x=d[,1],y=d[,3],xout = 0.1)$y)
      })      
})
saveRDS(realfdr,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/mast/plotdata/realfdr.rds')
saveRDS(reportfdr,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/mast/plotdata/reportfdr.rds')
