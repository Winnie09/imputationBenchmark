setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
allmtd = read.table('./resource/impute_method.txt',as.is = T)
allmtd = allmtd[,1]
mtd = 'scVI'
cl = 'H1975'
res = readRDS(file=paste0('./result/perf/sc_10x_5cl/',mtd,'.rds'))

sexpr = readRDS(paste0('./result/procimpute/cellbench/',mtd,'/sc_10x_5cl.rds'))
sexpr2 = readRDS(paste0('./result/procimpute/GSE81861/',mtd,'/GSE81861_Cell_Line_COUNT.rds'))

tab = table(c(paste0('10x_',sub('.*:','', colnames(sexpr))), paste0('fluidigm_',sub('_.*','',colnames(sexpr2)))))

bexpr_10x = readRDS('./data/bulkrna/cellbench/GSE86337_processed_count_average_replicates.rds')

raw = readRDS('./result/procimpute/cellbench/raw/sc_10x_5cl.rds')
rawcl = sub('.*:','',colnames(raw))


pseudobulk = rowMeans(raw[,rawcl == cl])
intergene = intersect(names(pseudobulk),rownames(bexpr_10x))
pd1 = data.frame(sc = pseudobulk[intergene], bulk = bexpr_10x[intergene,cl])

rawcor = cor(pseudobulk[intergene], bexpr_10x[intergene,cl], method='spearman')
intergene2 = intersect(rownames(sexpr), rownames(bexpr_10x))
pd2 = data.frame(sc = sexpr[intergene2, which(sub('.*:','', colnames(sexpr)) == cl)[1]], bulk = bexpr_10x[intergene,cl])

hmdf1 <- sapply(allmtd, function(mtd){
  res = readRDS(file=paste0('./result/perf/hm_cellline_cor/',mtd,'.rds'))
  sapply(res, median)
})
rownames(hmdf1) = paste0('fluidigm_',rownames(hmdf1))
hmdf2 <- sapply(allmtd, function(mtd){
  print(mtd)
  if (file.exists(paste0('./result/perf/sc_10x_5cl/',mtd,'.rds'))){
    res = readRDS(file=paste0('./result/perf/sc_10x_5cl/',mtd,'.rds'))
    sapply(res, median)      
  } 
})
rownames(hmdf2) = paste0('10x_', rownames(hmdf2))

# hmdf3 <- sapply(allmtd, function(mtd){ ####
#   print(mtd)
#   if (file.exists(paste0('./result/perf/hca/',mtd,'.rds'))){
#     res = readRDS(file=paste0('./result/perf/hca/',mtd,'.rds'))
#     sapply(res, median, na.rm=T)      
#   } 
# })
# hmdf3 = t(do.call(rbind, hmdf3))  ####


v_fluidigm = colMeans(hmdf1)
v_10x = colMeans(hmdf2)
# v_hca = colMeans(hmdf3) ####
saveRDS(v_fluidigm, './result/perf/assess/imp_eval_fluidigm.rds')
saveRDS(v_10x, './result/perf/assess/imp_eval_10x.rds')
# saveRDS(v_hca, './result/perf/assess/imp_eval_hca.rds') ####

hmdf = cbind(t(hmdf1),t(hmdf2))
library(reshape2)
hmdf = melt(hmdf) 
# tmpdf = melt(t(hmdf3)) ####
# hmdf = rbind(hmdf, tmpdf) ####

colnames(hmdf) <- c('method','ct','cor')
stat = tapply(hmdf$cor, hmdf$method, mean)
mtdorder = names(sort(tapply(hmdf$cor, hmdf$method, mean)))
saveRDS(mtdorder,'./result/perf/rank/imp_eval.rds')
saveRDS(stat, './result/perf/assess/imp_eval.rds')

hmdf$method = factor(hmdf$method, levels = mtdorder)
hmdf$ct = factor(as.character(hmdf$ct), levels=as.character(rev(sort(unique(hmdf$ct)))))
hmdf$ct = sapply(as.character(hmdf$ct), function(i) paste0(i,'(',tab[match(i,names(tab))],')') )
pd4 <- hmdf

df <- sapply(allmtd, function(mtd){
  res = readRDS(file=paste0('./result/perf/sc_10x_5cl/',mtd,'.rds'))
  res[[cl]]
})

pd = melt(df)
colnames(pd) = c('sc','method','cor')
pd$method = factor(as.character(pd$method), levels = mtdorder)
pd3 = pd
saveRDS(list(pd1=pd1,pd2=pd2,pd3=pd3,pd4=pd4),'./plot/pd/11_imp_eval_4subfig.pd.rds')
