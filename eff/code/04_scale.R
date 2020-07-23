setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
tb = readRDS('./eff/result/time.rds')
t = matrix(as.numeric(tb[,-1]),nrow=nrow(tb))
dimnames(t) = list(tb[,1],colnames(tb)[2:5])
t.bak = t
t[is.na(t)] = 60*24*3 + 1
score = matrix(0,nrow=nrow(t),ncol=ncol(t))
dimnames(score) = dimnames(t)
for (i in 1:ncol(t)){
  max = max(t[,i],na.rm=T)
  min = min(t[,i],na.rm=T)
  score[,i] = 1- (t[,i]-min)/(max-min)
}
saveRDS(score,'./eff/result/time_scaled.rds')
saveRDS(rowMeans(score),'./result/perf/assess/time.rds')

library(reshape2)
t = t.bak
colnames(t) <- c(1000,5000,50000,100000)
timedata <- melt(t)
colnames(timedata) <- c('method','cellnum','time')
timedata$cellnum <- log10(timedata$cellnum)
v <- sapply(unique(timedata$method),function(i){
  summary(lm(time~cellnum,data=timedata[timedata$method==i,]))$coefficients[2,1]  
})
names(v) = unique(timedata$method)
saveRDS(v,'./result/perf/assess/scalability.rds')

tb = readRDS('./eff/result/memory.rds')
tb = sub('K','',tb)
t = matrix(as.numeric(tb[,-1]),nrow=nrow(tb))
dimnames(t) = list(tb[,1],colnames(tb)[2:5])
t[is.na(t)] = max(t,na.rm = T)
score = matrix(0,nrow=nrow(t),ncol=ncol(t))
dimnames(score) = dimnames(t)
for (i in 1:ncol(t)){
  max = max(t[,i],na.rm=T)
  min = min(t[,i],na.rm=T)
  score[,i] = 1- (t[,i]-min)/(max-min)
}
saveRDS(score,'./eff/result/memory_scaled.rds')
saveRDS(rowMeans(score),'./result/perf/assess/memory.rds')

