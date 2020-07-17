allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/')
getf <- sub('.rds','',list.files("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/raw/"))
runf <- setdiff(allf, getf)
for (f in allf){
  raw = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',f, '/genebycell.rds'))
  saveRDS(raw,file=paste0("/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/raw/",f,'.rds'))  
}

