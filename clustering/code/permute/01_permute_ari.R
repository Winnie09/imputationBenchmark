source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
suppressMessages(library(mclust))
allf = c("sc_10x","sc_celseq2","sc_dropseq","sc_10x_5cl", "sc_celseq2_5cl_p1", "sc_celseq2_5cl_p2","sc_celseq2_5cl_p3")
f = allf[1]
set.seed(12345)
res <- sapply(allf,function(f){
  d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/raw/',f,'.rds'))
  ct = sub('.*:','',colnames(d))
  table(ct)
  ari = NULL
  for (s in seq(1,1e4)){
    set.seed(s)
    clu = get_permutation_clu(num.observation = length(ct), num.clu = length(unique(ct)))
    ari = c(ari,adjustedRandIndex(ct, clu))
  }
  ari
})
rm(d)
### pbmc
d = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/raw/sorted.rds')
ct = sub(':.*','',colnames(d))
rm(d)
ari = NULL
for (s in seq(1,1e4)){
  set.seed(s)
  clu = get_permutation_clu(num.observation = length(ct), num.clu = length(unique(ct)))
  ari = c(ari,adjustedRandIndex(ct, clu))
}

Res = cbind(res, ari)
colnames(Res) = c(allf, 'pbmc')
saveRDS(Res, '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/permute/permute_ari.rds')
