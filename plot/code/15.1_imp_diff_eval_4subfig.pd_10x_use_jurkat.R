setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
f1 <- list.files("./diff/result/10xcellline/")
f2 <- list.files("./diff/result/hm_cl/")
f <- intersect(f1,f2)
f <- sub('.rds','',f)
res <- NULL
hmdf <- NULL
for (sf in f) {
  set.seed(12345)
  d1 <- readRDS(paste0("./diff/result/10xcellline/",sf,'.rds'))
  d2 <- readRDS(paste0("./diff/result/hm_cl/",sf,'.rds'))
  d1 = list(as.vector(d1))
  names(d1) <- '293T_jurkat'
  names(d2) <- names(d2)
  d <- c(d1,d2)
  
  hmdf = rbind(hmdf, data.frame(Method=sf,Dataset=names(d),corMedian=sapply(d,median,na.rm=T),stringsAsFactors = F))
  d = lapply(d,sample,1e2)
  for (ds in names(d))
    res <- rbind(res,data.frame(Method=sf,Dataset=ds,Correlation=d[[ds]],stringsAsFactors = F))
}

platform = ifelse(hmdf[,'Dataset']=="293T_jurkat",'10x','fluidigm')
tmp = hmdf[platform=='10x', ]
v_10x = tapply(tmp$corMedian,list(tmp$Method),mean,na.rm=T)
tmp = hmdf[platform=='fluidigm', ]
v_fluidigm = tapply(tmp$corMedian,list(tmp$Method),mean,na.rm=T)
saveRDS(v_10x, './result/perf/assess/imp_diff_eval_10x.rds')
saveRDS(v_fluidigm, './result/perf/assess/imp_diff_eval_fluidigm.rds')
stat = rowMeans(tapply(hmdf[,'corMedian'], list(hmdf[,'Method'], platform), median, na.rm=T))
stat= stat[!is.na(stat)]
mtdorder = names(sort(stat))
hmdf = hmdf[complete.cases(hmdf),]
##
saveRDS(mtdorder,'./result/perf/rank/imp_diff_eval.rds')
saveRDS(stat,'./result/perf/assess/imp_diff_eval.rds')

library(reshape2)
mtd = 'saver'
sexpr = readRDS(paste0('./result/procimpute/10xcellline/',mtd,'/hg19.rds'))
sexpr_ct = sub('_.*','',colnames(sexpr))
###
bexpr_10x = readRDS('./data/bulkrna/expr/jurkat_hex.rds')
v1 = rowMeans(bexpr_10x[,which(colnames(bexpr_10x)=='239T')])
v2 = rowMeans(bexpr_10x[,which(colnames(bexpr_10x)=='Jurkat')])
bexpr_10x = cbind(v1,v2)
colnames(bexpr_10x) = c('293T','jurkat')
raw = readRDS('./result/procimpute/10xcellline/raw/hg19.rds')
rawcl = sub('_.*','',colnames(raw))
intergene = intersect(rownames(raw),rownames(bexpr_10x))
pbdiff = rowMeans(raw[intergene,rawcl == 'jurkat']) - rowMeans(raw[intergene,rawcl=='293T'])
pd1 <- data.frame(sc = pbdiff, bulk = bexpr_10x[intergene,2] - bexpr_10x[intergene,1])

rawcor = cor(pbdiff, (bexpr_10x[intergene,2] - bexpr_10x[intergene,1]), method='spearman')
intergene2 = intersect(rownames(sexpr), rownames(bexpr_10x))
impdiff = sexpr[intergene2,which(sexpr_ct=='jurkat')[1]] - sexpr[intergene2,which(sexpr_ct=='293T')[1]]
pd2 <- data.frame(sc = impdiff, bulk = bexpr_10x[intergene2,2] - bexpr_10x[intergene2,1])

pd = res[res[,'Dataset']=='293T_jurkat',]
pd[,'Method'] = factor(as.character(pd[,'Method']),levels=mtdorder)
pd3 <- pd

hmdf[,'Method'] = factor(hmdf[,'Method'], levels = mtdorder)
pd4 <- hmdf
saveRDS(list(pd1=pd1,pd2=pd2,pd3=pd3,pd4=pd4),'./plot/pd/15_imp_diff_eval_4subfig.pd.rds')

