mtd = commandArgs(trailingOnly = T)[1]
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/GSE81861/diff/wilcox/', mtd,'/res/')
dir.create(rdir, showWarnings = F, recursive = T)
imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/',mtd,'/GSE81861_Cell_Line_COUNT.rds'))
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/genebycell.rds')
imp = imp[,colnames(raw)]
g = intersect(rownames(imp),rownames(raw))
ct = sub('_.*','',colnames(raw))
imp = imp[g,ct=='GM12878']
raw = raw[g,ct=='GM12878']

df = expand.grid(c(10,20,30,40),c(10,20,30,40))
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

