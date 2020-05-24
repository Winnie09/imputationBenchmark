setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
f1 <- list.files("./diff/result/10xcellline/")
f2 <- list.files("./diff/result/hm_cl/")
f <- intersect(f1,f2)
f <- sub('.rds','',f)
res <- NULL

mtd = 'scVI'
sexpr = readRDS(paste0('./result/procimpute/cellbench/',mtd,'/sc_10x_5cl.rds'))
cl_10x <- sub('.*:','',colnames(sexpr))
n_10x <- table(cl_10x)

sexpr_fd <- readRDS(paste0('./result/procimpute/GSE81861/',mtd,'/GSE81861_Cell_Line_COUNT.rds'))
cl_fd <- sub('_.*','', colnames(sexpr_fd))
n_fd <- table(cl_fd)


hmdf <- NULL
for (sf in f) {
  set.seed(12345)
  d1 <- readRDS(paste0("./diff/result/sc_10x_5cl/",sf,'.rds'))
  d2 <- readRDS(paste0("./diff/result/hm_cl/",sf,'.rds'))
  d <- c(d1,d2)
  d = lapply(d,sample,1e2)
  tmpdf = rbind(data.frame(Method=sf,Dataset=names(d1),corMedian=sapply(d1,median,na.rm=T), platform='10x',stringsAsFactors = F),data.frame(Method=sf,Dataset=names(d2),corMedian=sapply(d2,median,na.rm=T), platform='fluidigm',stringsAsFactors = F))
  hmdf = rbind(hmdf, tmpdf)
  

   for (ds in names(d))
    res <- rbind(res,data.frame(Method=sf,Dataset=ds,Correlation=d[[ds]],platform=ifelse(ds%in%names(d1),'10x','fluidigm'),stringsAsFactors = F))
}
hmdf = hmdf[complete.cases(hmdf),]
nm_10x <- NULL
for (i in unique(hmdf[hmdf[,'platform']=='10x','Dataset'])){
  tmp1 <- sub('_.*','',i)
  tmp2 <- sub('.*_','',i)
  tmp <- paste0(tmp1, '(', n_10x[tmp1], ')_',tmp2, '(', n_10x[tmp2], ')')
  hmdf[hmdf[,'platform']=='10x' & hmdf[,'Dataset']==i, 'Dataset'] <-tmp
  nm_10x <- c(nm_10x,mean(c(n_10x[tmp1], n_10x[tmp2])))
  names(nm_10x)[length(nm_10x)] <- tmp
}
            
nm_fd <- NULL
for (i in unique(hmdf[hmdf[,'platform']=='fluidigm','Dataset'])){
  tmp1 <- sub('_.*','',i)
  tmp2 <- sub('.*_','',i)
  tmp <- paste0(tmp1, '(', n_fd[tmp1], ')_',tmp2, '(', n_fd[tmp2], ')')
  hmdf[hmdf[,'platform']=='fluidigm' & hmdf[,'Dataset']==i, 'Dataset'] <-tmp
    nm_fd <- c(nm_fd,mean(c(n_fd[tmp1], n_fd[tmp2])))
  names(nm_fd)[length(nm_fd)] <- tmp
}

o1 <- order(nm_10x, decreasing=TRUE)
o2 <- order(nm_fd, decreasing=TRUE)
o <- c(o1, o2 + length(o1))
hmdf[,'Method'] = factor(hmdf[,'Method'], levels = mtdorder)
hmdf[,'platform'] = factor(hmdf[,'platform'], levels = c('10x','fluidigm'))
hmdf[,'Dataset'] = factor(hmdf[,'Dataset'], levels = c(names(nm_10x), names(nm_fd))[o])
pd4 <- hmdf

tmp = hmdf[hmdf$platform=='10x',]
v_10x = tapply(tmp$corMedian,list(tmp$Method),mean,na.rm=T)
tmp = hmdf[hmdf$platform=='fluidigm', ]
v_fluidigm = tapply(tmp$corMedian,list(tmp$Method),mean,na.rm=T)
saveRDS(v_10x, './result/perf/assess/imp_diff_eval_10x.rds')
saveRDS(v_fluidigm, './result/perf/assess/imp_diff_eval_fluidigm.rds')
stat = rowMeans(tapply(hmdf[,'corMedian'], list(hmdf[,'Method'], hmdf$platform), median, na.rm=T))
stat= stat[!is.na(stat)]
mtdorder = names(sort(stat))

saveRDS(mtdorder,'./result/perf/rank/imp_diff_eval.rds')
saveRDS(stat,'./result/perf/assess/imp_diff_eval.rds')

library(reshape2)
mtd = 'scVI'
sexpr = readRDS(paste0('./result/procimpute/cellbench/',mtd,'/sc_10x_5cl.rds'))
sexpr_ct = sub('.*:','',colnames(sexpr))
###
bexpr_10x = readRDS('./data/bulkrna/cellbench/GSE86337_processed_count_average_replicates.rds')
bexpr_10x = bexpr_10x[,c('H1975','H2228')]

raw = readRDS('./result/procimpute/cellbench/raw/sc_10x_5cl.rds')
rawcl = sub('.*:','',colnames(raw))
intergene = intersect(rownames(raw),rownames(bexpr_10x))
pbdiff = rowMeans(raw[intergene,rawcl == 'H1975']) - rowMeans(raw[intergene,rawcl=='H2228'])
pd1 <- data.frame(sc = pbdiff, bulk = bexpr_10x[intergene,1] - bexpr_10x[intergene,2])

rawcor = cor(pbdiff, (bexpr_10x[intergene,1] - bexpr_10x[intergene,2]), method='spearman')
intergene2 = intersect(rownames(sexpr), rownames(bexpr_10x))
impdiff = sexpr[intergene2,which(sexpr_ct=='H1975')[1]] - sexpr[intergene2,which(sexpr_ct=='H2228')[1]]
pd2 <- data.frame(sc = impdiff, bulk = bexpr_10x[intergene2,1] - bexpr_10x[intergene2,2])

pd = res[res[,'Dataset']=='H2228_H1975',]
pd[,'Method'] = factor(as.character(pd[,'Method']),levels=mtdorder)
pd3 <- pd

saveRDS(list(pd1=pd1,pd2=pd2,pd3=pd3,pd4=pd4),'./plot/pd/15_imp_diff_eval_4subfig.pd.rds')
