setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
allmtd = read.table('./resource/impute_method.txt',as.is = T)
allmtd = allmtd[,1]
#allmtd = setdiff(allmtd,'viper')

mtd = 'saver'
cl = '293T'
res = readRDS(file=paste0('./result/perf/hm_cellline_cor/',mtd,'.rds'))
sexpr = readRDS(paste0('./result/procimpute/10xcellline/',mtd,'/hg19.rds'))
sexpr2 = readRDS(paste0('./result/procimpute/GSE81861/',mtd,'/GSE81861_Cell_Line_COUNT.rds'))
tab = table(sub('_.*','',c(colnames(sexpr), colnames(sexpr2))))
##c = sum(grepl('293T',colnames(sexpr))), sum(grepl('jurkat',colnames(sexpr)))
# x: average all raw single cells in cl, 
# y: bulk cl. each dot is a gene
bexpr_10x = readRDS('./data/bulkrna/expr/jurkat_hex.rds')
v1 = rowMeans(bexpr_10x[,which(colnames(bexpr_10x)=='239T')])
v2 = rowMeans(bexpr_10x[,which(colnames(bexpr_10x)=='Jurkat')])
bexpr_10x = cbind(v1,v2)
colnames(bexpr_10x) = c('293T','jurkat')

raw = readRDS('./result/procimpute/10xcellline/raw/hg19.rds')
rawcl = sub('_.*','',colnames(raw))
pseudobulk = rowMeans(raw[,rawcl == cl])
intergene = intersect(names(pseudobulk),rownames(bexpr_10x))
pd1 = data.frame(sc = pseudobulk[intergene], bulk = bexpr_10x[intergene,cl])


rawcor = cor(pseudobulk[intergene], bexpr_10x[intergene,cl], method='spearman')
intergene2 = intersect(rownames(sexpr), rownames(bexpr_10x))
pd2 = data.frame(sc = sexpr[intergene2, which(sub('_.*','', colnames(sexpr)) == cl)[1]], bulk = bexpr_10x[intergene,cl])

hmdf1 <- sapply(allmtd, function(mtd){
  res = readRDS(file=paste0('./result/perf/hm_cellline_cor/',mtd,'.rds'))
  sapply(res, median)
})
hmdf2 <- sapply(allmtd, function(mtd){
  print(mtd)
  if (file.exists(paste0('./result/perf/10xcellline_cor/',mtd,'.rds'))){
    res = readRDS(file=paste0('./result/perf/10xcellline_cor/',mtd,'.rds'))
    sapply(res, median)      
  } 
})
v_fluidigm = colMeans(hmdf1)
v_10x = colMeans(hmdf2)
saveRDS(v_fluidigm, './result/perf/assess/imp_eval_fluidigm.rds')
saveRDS(v_10x, './result/perf/assess/imp_eval_10x.rds')

hmdf = cbind(t(hmdf1),t(hmdf2))
library(reshape2)
hmdf = melt(hmdf)
colnames(hmdf) <- c('method','ct','cor')
stat = tapply(hmdf$cor, hmdf$method, mean)
mtdorder = names(sort(tapply(hmdf$cor, hmdf$method, mean)))
saveRDS(mtdorder,'./result/perf/rank/imp_eval.rds')
saveRDS(stat, './result/perf/assess/imp_eval.rds')

hmdf$method = factor(hmdf$method, levels = mtdorder)
hmdf$ct = factor(as.character(hmdf$ct), levels=c("293T","jurkat","A549","GM12878","H1","IMR90","K562" ) )
hmdf$ct = sapply(as.character(hmdf$ct), function(i) paste0(i,'(',tab[match(i,names(tab))],')') )
pd4 <- hmdf

df <- sapply(allmtd, function(mtd){
  res = readRDS(file=paste0('./result/perf/10xcellline_cor/',mtd,'.rds'))
  res[[cl]]
})

pd = melt(df)
colnames(pd) = c('sc','method','cor')
pd$method = factor(as.character(pd$method), levels = mtdorder)
pd3 = pd
saveRDS(list(pd1=pd1,pd2=pd2,pd3=pd3,pd4=pd4),'./plot/pd/11_imp_eval_4subfig.pd.rds')
