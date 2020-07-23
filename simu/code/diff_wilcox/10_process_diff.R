mtd = commandArgs(trailingOnly = T)[1]
allf = sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_wilcox/diff/',mtd,'/')))
res = list()
for (f in allf){
      print(f)
      res[[f]] = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_wilcox/diff/',mtd,'/',f,'.rds'))
}
saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_wilcox/procdiff/',mtd,'.rds'))
