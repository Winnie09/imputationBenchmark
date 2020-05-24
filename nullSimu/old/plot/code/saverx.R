mtd = 'saverx'
print(mtd)
allf =  list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/',mtd))
a = sapply(allf, function(f){
  mat = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/',mtd,'/', f))
  library(ggplot2)
  mat = mat[rowSums(mat)>0,]
  res <- prcomp(t(mat), scale=T)$x
  v = apply(res,2,var)
  dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',mtd), showWarnings = F)
  saveRDS(list(pc = res[,1:10], var = v),paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',mtd,'/',f))
  return(0)      
})

