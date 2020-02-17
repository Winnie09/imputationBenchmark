mtd = commandArgs(trailingOnly = T)[1]
suppressMessages(library(MAST))
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/sc_10x_5cl/diff/mast/', mtd,'/res/')
dir.create(rdir, showWarnings = F, recursive = T)
imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/sc_10x_5cl.rds'))
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/sc_10x_5cl/genebycell.rds')
imp = imp[,colnames(raw)]
g = intersect(rownames(imp),rownames(raw))
ct = sub('.*:','',colnames(raw))
imp = imp[g,ct=='A549']
raw = raw[g,ct=='A549']

df = expand.grid(c(10,50,100,500),c(10,50,100,500))
colnames(df) =  c('n1','n2')
df = df[df[,'n1']<=df[,'n2'],]
###
for (i in 1:nrow(df)){
  cn1 = df[i,'n1']
  cn2 = df[i,'n2']
  set.seed(12345)
  id = sample(1:ncol(imp), cn1+cn2)
  id1 = id[1:cn1]
  id2 = setdiff(id,id1)
  expr = imp[,id]
  count = raw[,id]
  cdr <- scale(colMeans(count > 0))
  cluster <- rep('clu1',ncol(expr)) ## length = #cells
  cluster[(cn1+1):ncol(expr)] <- 'clu2' ###
  sca <- FromMatrix(expr,data.frame(wellKey=colnames(expr),cluster=cluster,cngeneson=cdr), data.frame(primerid=row.names(expr),Gene=row.names(expr)), check_sanity = FALSE)
  zlmCond <- zlm(~cluster+cngeneson, sca)
  summaryCond <- summary(zlmCond, doLRT="clusterclu2")
  summaryDt <- as.data.frame(summaryCond$datatable)
  pval <- summaryDt[summaryDt$component == "H" & summaryDt$contrast == "clusterclu2",c(1,4)]
  lfc <- summaryDt[summaryDt$component == "logFC" & summaryDt$contrast == "clusterclu2",c(1,7)]
  combine <- merge(pval,lfc)
  colnames(combine) <- c("Gene","pvalue","Log-foldchange")
  row.names(combine) <- combine[,1]
  res <- combine
  colnames(res)[colnames(res)=='Log-foldchange'] <- "stat"
  res <- res[order(res[,'pvalue'],-abs(res[,'stat'])),]
  saveRDS(res, paste0(rdir,cn1,'_',cn2,'.rds'))  
}


