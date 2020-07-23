DEMethod = commandArgs(trailingOnly = T)[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_', DEMethod,'/diff'))
allf = setdiff(allf,'zinbwave')
d = list()
for (i in sub('.rds', '', allf)){
  file = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',i,'.rds')
  if (file.exists(file)){
    tmp = readRDS(file)
    d[[i]] = tmp  
  }
  
}
res <- NULL
for (i in names(d)) {
  print(i)
  for (j in names(d[[i]])) {
    tmp <- d[[i]][[j]]
    bound <- approx(x=tmp[,2],y=tmp[,1],xout=0.25)$y
    tmp <- rbind(tmp[tmp[,2] < 0.25,1:2],c(bound,0.25))
    area <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm=T)/0.25
    res <- rbind(res,data.frame(i,sub('.rds','',j),area,stringsAsFactors = F))
  }
}
colnames(res) <- c('method','setting','AUC')
allres <- res
res$cellnum <- paste0(strsplit(sub('.rds','',res[,2]), '_')[[1]][1],'_',strsplit(sub('.rds','',res[,2]), '_')[[1]][2])
res[,1] <- factor(res[,1],levels=names(sort(tapply(res$AUC,res$method,median))))
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plotdata/'), showWarnings = F)
saveRDS(res, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plotdata/auc.rds'))
