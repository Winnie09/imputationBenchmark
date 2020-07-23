allf =  list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/')
for (f in allf){
      raw = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/',f,'/norm_genebycell.rds'))
      saveRDS(raw,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/raw/',f,'.rds'))
}

