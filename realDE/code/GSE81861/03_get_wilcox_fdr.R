mtd = commandArgs(trailingOnly = T)[1]
suppressMessages(library(MAST))
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/diff/wilcox/', mtd,'/res/')
dir.create(rdir, showWarnings = F, recursive = T)
imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/',mtd,'/GSE81861_Cell_Line_COUNT.rds'))
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/genebycell.rds')
raw = raw[,colnames(imp)]
g = intersect(rownames(imp), rownames(raw))
imp = imp[g, ]
raw = raw[g, ]
ct = sub('_.*','',colnames(raw))
allf = sub('_diffgene.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/GSE81861/bulkdiff/'))

for (f in allf){
  ct1 = sub('_.*','',f)
  ct2 = sub('.*_','',f)
  pval <- sapply(1:nrow(imp), function(i) {
    v1 = imp[i, ct == ct1]
    v2 = imp[i, ct == ct2]
    if (length(unique(c(v1,v2)))==1){
      1
    } else {
      wilcox.test(v1, v2)$p.value  
    }
  })
  fdrv <- p.adjust(pval,method='fdr')
  df = data.frame(Gene=rownames(imp), fdr = fdrv)
  saveRDS(df, paste0(rdir,f,'.rds'))  
}

