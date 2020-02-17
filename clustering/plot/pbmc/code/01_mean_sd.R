allmtd = readLines('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/resource/impute_method_latent.txt')
for (mtd in allmtd){
  print(mtd)
  if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/',mtd,'/sorted.rds'))){
    d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/',mtd,'/sorted.rds'))
    library(splines)
    m = rowMeans(d)
    v = apply(d,1,sd)
    pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/plot/mean_sd_originalScale/',mtd,'.pdf'),width=3,height=3)
    plot(v~m, pch=20,cex=.5,main=paste(mtd,round(cor(m,v),2)), xlab='mean',ylab='sd')
    dev.off()
    ## log scale
    m = log2(m)
    v = log2(v)
    pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/plot/mean_sd/',mtd,'.pdf'),width=3,height=3)
    plot(v~m, pch=20,cex=.5,main=paste(mtd,round(cor(m,v),2)), xlab='log2(mean)',ylab='log2(sd)')
    dev.off()
  }
}

