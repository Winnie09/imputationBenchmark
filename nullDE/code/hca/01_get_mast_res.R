mtch = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cell_celltype_cluster.rds')
cluct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cluster_celltype.rds')
cluct = cluct[complete.cases(cluct),]
cluct_s = paste0(cluct[,1],'_',cluct[,2])
mtch_s = paste0(mtch$cluster, '_', mtch$celltype)
select_cluct = names(sort(table(mtch_s[mtch_s%in%cluct_s]),decreasing=T)[1])
select_cell = mtch[mtch$cluster == sub('_.*','',select_cluct) & mtch$celltype == sub('.*_','',select_cluct), 'cell']


mtd = as.character(commandArgs(trailingOnly = T)[[1]])
suppressMessages(library(MAST))
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/hca/diff/mast/', mtd,'/res/')
dir.create(rdir, showWarnings = F, recursive = T)
imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd,'/MantonBM6.rds'))
if (grepl('A.1',colnames(imp)[1])){
  colnames(imp) = sub('\\.','-',colnames(imp))
}
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/MantonBM6/genebycell.rds')
g = intersect(rownames(imp), rownames(raw))
imp = imp[g,select_cell]
raw = raw[g,select_cell]

df = expand.grid(c(30,50,70,90),c(30,50,70,90))
colnames(df) =  c('n1','n2')
df = df[df[,'n1']<=df[,'n2'],]

for (i in 1:nrow(df)){
  cn1 = df[i,'n1']
  cn2 = df[i,'n2']
  set.seed(12345)
  id = sample(1:ncol(imp), cn1+cn2)
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

