mtd = commandArgs(trailingOnly = T)[1]
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/sc_10x_5cl/diff/wilcox/', mtd,'/res/')
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

for (i in 1:nrow(df)){
  print(i)
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
  
  wval <- sapply(1:nrow(expr), function(i) {
    v1 = expr[i, 1:cn1]
    v2 = expr[i, ((cn1+1):ncol(expr))]
    if (length(unique(c(v1,v2)))==1){
      1
    } else {
      if (mean(v1) >= mean(v2)){
        wilcox.test(v1, v2)$statistic  
      } else {
        wilcox.test(v2, v1)$statistic  
      }
      
    }
  })
  
  fdrv <- p.adjust(pval,method='fdr')
  res = data.frame(Gene=rownames(expr), fdr = fdrv, statistics=wval, pval = pval)
  res = res[order(res[,'fdr'], -res[,'statistics']),]
  print(summary(res[,2]))
  saveRDS(res, paste0(rdir,cn1,'_',cn2,'.rds'))
}


