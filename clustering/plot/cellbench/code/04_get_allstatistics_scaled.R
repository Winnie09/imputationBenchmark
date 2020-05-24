clumethod = as.character(commandArgs(trailingOnly = T)[1])
pm = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/permute/permute_ari.rds')
min.ari = apply(pm,2,min)[1:7]

res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/',clumethod,'_cellbench_allstatistics.rds'))
nclu = as.numeric(sapply(res$data,function(i) ifelse(grepl('5',i),'5','3')))
min.ari.tmp = min.ari[as.character(res$data)]

mat = sapply(colnames(res)[3:6],function(i){
  if (i == 'Hacc' || i == 'Hpur'){
    max = log(nclu)
    min = 0
    v = 1-(res[,i] - min)/(max - min)
  } else if (i == 'ARI'){
    max = 1
    min = min.ari.tmp
    v = (res[,i] - min.ari.tmp)/(max - min.ari.tmp)
  } else if (i == 'medianSil'){
    max = 1
    min = -1
    v = (res[,i] - min)/(max - min)
  }
})
mat = data.frame(method=as.character(res[,1]),data=as.character(res[,2]),mat)
saveRDS(mat,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/',clumethod,'_cellbench_allstatistics_scaled.rds'))


mat2 = sapply(colnames(mat)[3:6],function(i){
  tapply(mat[,i],list(res[,1]),mean,na.rm=T)
})
saveRDS(mat2,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/stat/',clumethod,'_cellbench_statistics_scaled_summary.rds'))

