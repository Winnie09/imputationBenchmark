library(mclust)
library(cluster)
# allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/')
# allmtd = setdiff(allmtd,c('deepimpute','scimpute'))
mtd = as.character(commandArgs(trailingOnly = T)[1])
print(mtd)
allf = c("sc_10x","sc_celseq2","sc_dropseq","sc_10x_5cl", "sc_celseq2_5cl_p1", "sc_celseq2_5cl_p2","sc_celseq2_5cl_p3")
df = data.frame(allf = allf, k = c(rep(3,3), rep(5,4)))
existf = sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/',mtd,'/')))
if (length(existf)>0){
  ACC <- PUR <- ARI <- meanSil <- NULL
  for (f in intersect(allf, existf)){
    print(f)
    clu = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/',mtd,'/',f,'.rds'))
    ct = sub('.*:','',names(clu))
    acc <- -mean(sapply(unique(clu),function(i){
      p = table(ct[clu==i])/ sum(clu==i)
      sum(p * log(p))
    }))
    pur <- -mean(sapply(unique(ct),function(sct){
      p = table(clu[ct==sct])/ sum(ct == sct)
      sum(p * log(p))
    }))
    ACC = c(ACC, acc)
    PUR = c(PUR, pur)  
    ## adjusted rank index
    suppressMessages(library(mclust))
    ARI <- c(ARI,adjustedRandIndex(ct,clu))
    ## sum of silhouette
    mat = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/',f,'.rds'))
    d = dist(t(mat))    ## use cell by gene for dist()
    s = silhouette(clu, dist=d)
    meanSil = c(meanSil,mean(s[,3]))
  }
  df = data.frame(method=mtd, data = intersect(allf, existf), Hacc =ACC, Hpur=PUR, ARI = ARI, meanSil = meanSil)
}
saveRDS(df,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/meanSil/',mtd,'.rds'))

