mtch = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cell_celltype_cluster.rds')
cluct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cluster_celltype.rds')
cluct = cluct[complete.cases(cluct),]
cluct_s = paste0(cluct[,1],'_',cluct[,2])
mtch_s = paste0(mtch$cluster, '_', mtch$celltype)
select_cluct = names(sort(table(mtch_s[mtch_s%in%cluct_s]),decreasing=T)[1])
select_cell = mtch[mtch$cluster == sub('_.*','',select_cluct) & mtch$celltype == sub('.*_','',select_cluct), 'cell']


mtd = as.character(commandArgs(trailingOnly = T)[[1]])
suppressMessages(library(MAST))
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/hca/diff/wilcox/', mtd,'/res/')
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
  pval <- sapply(1:nrow(expr), function(i) {
    v1 = expr[i, 1:cn1]
    v2 = expr[i, ((cn1+1):ncol(expr))]
    if (length(unique(c(v1,v2)))==1){
      1
    } else {
      wilcox.test(v1, v2)$p.value  
    }
  })
  fdrv <- p.adjust(pval,method='fdr')
  res = data.frame(Gene=rownames(expr), fdr = fdrv)
  res = res[order(res[,'fdr']),]
  saveRDS(res, paste0(rdir,cn1,'_',cn2,'.rds'))  
}

