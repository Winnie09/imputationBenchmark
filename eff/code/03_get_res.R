timefunc <- function(i) {
  if (!grepl('-',i)) {
    t = as.numeric(times(i))
  } else {
    t = as.numeric(sub('-.*','',i))+as.numeric(times(sub('.*-','',i)))
  }
  round(60 * 24 * t,3)
}

setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/code/impute/')
af = list.files(getwd())
library(chron)
eff <- lapply(af,function(f){
  print(f)
  d=readLines(paste0(f,'/out'))
  a = sapply(d[grep('whou10',d)+1],function(i) {
    tmp <- strsplit(i,' ')[[1]]
    tmp <- tmp[nchar(tmp) > 0]
    if ('COMPLETED' %in% tmp) {
      tmp[(length(tmp)-1):length(tmp)]  
    } else if ('FAILED' %in% tmp | 'RUNNING' %in% tmp) {
      c(NA,NA)
    }
  },USE.NAMES = F)
  if (is.list(a)){
    a = a[!is.na(a)]
    a = do.call(cbind,a)
  }
  a = a[,order(as.numeric(sub('K','',a[2,])))]
  a[1,] <- sapply(a[1,],timefunc)
  a = cbind(f,a,matrix(NA,2,4-ncol(a)))
  colnames(a) = c('method','1k','5k','50k','100k')
  a
})
names(eff) = af
res = do.call(rbind,eff)


time = res[seq(1,nrow(res),2),]
mem = res[seq(2,nrow(res),2),]
saveRDS(time,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/result/time.rds')
saveRDS(mem,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/result/memory.rds')
