allmtd = readLines('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/resource/impute_method_latent.txt')
#getmtd = sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/pd/umap/'))
#runmtd = setdiff(allmtd,getmtd)
for (mtd in allmtd){
  print(mtd)
  af <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd))
  dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/pd/pc/', mtd), recursive = T)
  dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/pd/umap/',mtd), recursive = T)
  if (length(af)>=1){
    res <- sapply(af, function(f){
      print(f)
      if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/',f)) & !file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/pd/pc/',mtd,'/', f))){
        d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/',f))
        if (!grepl('latent', mtd)){
          library(splines)
          m = rowMeans(d)
          v = apply(d,1,sd)
          m = log2(m)
          v = log2(v)
          id <- intersect(which(abs(v)!=Inf), which(abs(m)!=Inf))
          resid = resid(lm(v[id]~bs(m[id],df=6)))
          d = d[names(sort(resid,decreasing = T))[1:1e3], ]
          d = prcomp(t(d),scale. = T)$x[,1:10] ## cell by PC
          
        } else {
          d <- t(d)
          
        }
        saveRDS(d,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/pd/pc/',mtd,'/', f))
        library(umap)
        set.seed(12345)
        umap=umap(d)$layout
        saveRDS(umap,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/pd/umap/',mtd,'/',f))
        rm(d)
        return(0)
      }
    })
  }
}

