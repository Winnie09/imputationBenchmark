rm(list=ls())
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
tb = read.csv('./doc/evalution_table.csv',stringsAsFactors = F)
allf = tb$evaluation
allmtd = readLines('./resource/impute_method_latent.txt')
mat <- matrix(NA,nrow=length(allmtd),ncol=length(allf),dimnames = list(allmtd,allf))
for (f in allf) {
  print(f)
  tmp = readRDS(paste0('./result/perf/assess/',f,'.rds'))
  if (f == 'scalability'){
    max = max(tmp,na.rm=T)
    min = min(tmp,na.rm=T)
    tmp = 1 - (tmp-min)/(max-min)
  }
  mat[names(tmp),f] <-tmp
}
  
mat[,grepl('nullDE', colnames(mat))] = 1- mat[,grepl('nullDE', colnames(mat))]/50500
tmp = mat[,1:14]
tmp[is.na(tmp)] <- 0
mat[,1:14] = tmp
cid = !grepl('trajectory',colnames(mat))
mat['saucie',cid] <- pmax(mat['saucie',cid],mat['saucie_latent',cid],na.rm = T)
mat['scscope',cid] <- pmax(mat['saucie',cid],mat['scscope_latent',cid],na.rm = T)
mat['scVI',cid] <- pmax(mat['scVI',cid],mat['scVI_latent',cid],na.rm = T)
mat = mat[!grepl('latent',rownames(mat)),]
saveRDS(mat, './plot/plot/all_eval_score_noNA.rds')

