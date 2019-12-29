allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute')
tmp <- sapply(allmtd,function(mtd){
  print(mtd)
  allf =  list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/',mtd))
  a = sapply(allf, function(f){
    print(f)
    mat = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/',mtd,'/', f))
    library(ggplot2)
    if (!grepl('latent',mtd)){
      mat = mat[rowSums(mat)>0,]
    }
    rsd <- apply(mat,1,sd)
    rm <- rowMeans(mat)
    cv <- rsd/rm
    mat <- mat[cv > median(cv),]
    res <- prcomp(t(mat), scale=T)$x
    v = apply(res,2,var)
    dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',mtd), showWarnings = F)
    saveRDS(list(pc = res[,1:min(10,ncol(res))], var = v),paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',mtd,'/',f))
    return(0)      
  })
  return(0)
})

