library(ggplot2)
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute')
len = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19.rds')
df = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg19_geneid_genename.rds')
names(len) = df[match(names(len), df$geneid), 'genename']
# allf =  sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/alra')))
# alld = list()
# for (f in allf){
#   alld[[f]] = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/',f,'/norm_genebycell.rds'))
# }
tmp <- sapply(allmtd,function(mtd){
  print(mtd)
  allf =  sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/',mtd)))
  a = sapply(allf, function(f){
    mat = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/',mtd,'/', f,'.rds'))
    mat = mat[rowSums(mat)>0, ]
    # d = alld[[f]]
    g = intersect(rownames(mat), names(len))
    mat = mat[g,]
    
    res <- prcomp(mat, scale=T)$x
    v = apply(res,2,var)
    geneExprMean = apply(mat,1,mean)
    geneExprVar = apply(mat,1,var)
    dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca_bygene/',mtd), showWarnings = F)
    saveRDS(list(pc = res[,1:10], pcvar = v, geneExprMean=geneExprMean, geneExprVar=geneExprVar, geneLen = len[g]),paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca_bygene/',mtd,'/',f,'.rds'))
    return(0)      
  })
  return(0)
})

