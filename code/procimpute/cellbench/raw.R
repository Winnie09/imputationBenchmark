allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/')
res <- lapply(allf, function(f){
      raw = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/', f,'/norm_genebycell.rds'))
      saveRDS(raw,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/raw/',f,'.rds'))
})
