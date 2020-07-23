allmtd = readLines('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/resource/impute_method_latent.txt')
#getmtd = sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/pd/umap/'))
#runmtd = setdiff(allmtd,getmtd)
for (mtd in allmtd){
  print(mtd)
  if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/',mtd,'/sorted.rds'))){
    d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/',mtd,'/sorted.rds'))
    if (grepl('latent', mtd)){
      d <- t(d)
    } else {
      library(splines)
      m = rowMeans(d)
      v = apply(d,1,sd)
      m = log2(m)
      v = log2(v)
      resid = resid(lm(v~bs(m,df=6)))
      d = d[names(sort(resid,decreasing = T))[1:1e3], ]
      d = prcomp(t(d),scale. = T)$x[,1:10] ## cell by PC
    }
    saveRDS(d,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/pd/pc/',mtd,'.rds'))
    library(umap)
    set.seed(12345)
    umap=umap(d)$layout
    saveRDS(umap,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/pd/umap/',mtd,'.rds'))
    rm(d)
  }
}

