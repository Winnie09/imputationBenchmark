
allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/check_cellline_differential/result/fdr/diffNumCell')
num = NULL
for (f in allf){
  tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/check_cellline_differential/result/fdr/diffNumCell/', f))
  num = c(num,tmp)
}
names(num) = allf
summary(num)


allf = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/check_cellline_differential/result/fdr/sameNumCell')
length(allf)
num = NULL
for (f in allf){
  tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/check_cellline_differential/result/fdr/sameNumCell/', f))
  num = c(num,tmp)
}
names(num) = allf
summary(num)
