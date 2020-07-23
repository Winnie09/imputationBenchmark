DEMethod = commandArgs(trailingOnly = T)[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/diff'))
allf = setdiff(allf,'zinbwave')
d = list()
for (i in sub('.rds', '', allf)){
  file = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',i,'.rds')
  if (file.exists(file)){
    tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',i,'.rds'))
    d[[i]] = tmp  
  }
}
res <- NULL
for (i in names(d)) {
  print(i)
  for (j in names(d[[i]])) {
    print(j)
    tmp <- d[[i]][[j]]
    if (length(unique(tmp[,3])) == 1){
      next
    }
    bound <- approx(x=tmp[,3],y=tmp[,2],xout=0.25)$y
    tmp <- rbind(tmp[tmp[,3] < 0.25,2:3],c(bound,0.25))
    tmp <- unique(tmp)
    # area <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm = T)-0.25*0.25/2
    area <- (0.25*0.25/2) - sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm = T)
    res <- rbind(res,data.frame(i,sub('.rds','',j),area,stringsAsFactors = F))
  }
}
colnames(res) <- c('method','setting','diff')
res[,1] <- factor(res[,1],levels=names(sort(tapply(res$diff,res$method,median))))
res$cellnum <- paste0(strsplit(sub('.rds','',res[,2]), '_')[[1]][1],'_',strsplit(sub('.rds','',res[,2]), '_')[[1]][2])

pddir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/', DEMethod, '/plotdata/')
dir.create(pddir, showWarnings = F)
saveRDS(res, paste0(pddir,'diff.rds'))

