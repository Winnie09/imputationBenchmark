mtd = as.character(commandArgs(trailingOnly = T)[1])
clumethod = as.character(commandArgs(trailingOnly = T)[2])
library(mclust)
library(cluster)
print(mtd)
print(clumethod)
existf = sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/pbmc/',clumethod,'/',mtd,'/')))
set.seed(12345)
if (length(existf)>0){
  ACC <- PUR <- ARI <- medianSil <- NULL
  for (f in existf){
    print(f)
    clu = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/pbmc/',clumethod,'/',mtd,'/',f,'.rds'))
    ct = sub(':.*','',names(clu))
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
    mat = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/',mtd,'/',f,'.rds'))
    sid = sample(1:ncol(mat), 10000)
    d = dist(t(mat[,sid]))    ## use cell by gene for dist()
    s = silhouette(clu[sid], dist=d)
    medianSil = c(medianSil,median(s[,3]))
  }
  df = data.frame(method=mtd, data = existf, Hacc =ACC, Hpur=PUR, ARI = ARI, medianSil = medianSil)
}
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/pbmc/',clumethod,'/medianSil/'),recursive=T,showWarnings = F)
saveRDS(df,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/perf/pbmc/',clumethod,'/medianSil/',mtd,'.rds'))
