DEMethod = commandArgs(trailingOnly = T)[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/diff'))
d = list()
for (i in sub('.rds', '', allf)){
  tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',i,'.rds'))
  d[[i]] = tmp
}
res <- NULL
for (i in names(d)) {
  print(i)
  for (j in names(d[[i]])) {
    tmp <- d[[i]][[j]]
    tmp = cbind(tmp, Precision = 1-tmp[,'Real_FDR'])
    #bound <- approx(y=tmp[,'Precision'],x=tmp[,'Sensitivity'],xout=0.25)$y
    #tmp <- rbind(tmp[tmp[,'Real_FDR'] < 0.25, c('Sensitivity','Real_FDR')],c(bound,0.25))
    area <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,'Precision']+tmp[i,'Precision'])*(tmp[i,'Sensitivity']-tmp[i-1,'Sensitivity'])/2),na.rm=T)
    res <- rbind(res,data.frame(i,sub('.rds','',j),area,stringsAsFactors = F))
  }
}
colnames(res) <- c('method','setting','AUC')
allres <- res
res$cellnum <- paste0(strsplit(sub('.rds','',res[,2]), '_')[[1]][1],'_',strsplit(sub('.rds','',res[,2]), '_')[[1]][2])
res[,1] <- factor(res[,1],levels=names(sort(tapply(res$AUC,res$method,median))))
saveRDS(res, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plotdata/precision_recall_auc.rds'))

