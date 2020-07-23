source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
allf <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/')
pattern = as.character(commandArgs(trailingOnly = T)[1])
print(pattern)
m = as.character(commandArgs(trailingOnly = T)[2])
allf = allf[grep(pattern,allf)]
cormat <- sapply(allf,function(f){
    print(f)
      dg = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f,'/diffgn.rds'))
      raw = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f,'/genebycell.rds'))
      dg = dg[dg%in%row.names(raw)]      
      if (m == 'raw'){
        sample(as.vector(corfunc(t(raw[dg,]), t(raw[setdiff(row.names(raw),dg),]))), 1e3)  
      } else {
        if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m,'/',f,'.rds'))) {
          d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m,'/',f,'.rds'))
          m1d = d[dg,]
          m1nd = d[setdiff(row.names(d),dg),]
          sample(as.vector(corfunc(t(m1d),t(m1nd))),1e3)
        }  
      }
  })
saveRDS(cormat, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/genecor/result/diff_ndiff_cor/', pattern,'_',m,'.rds'))


