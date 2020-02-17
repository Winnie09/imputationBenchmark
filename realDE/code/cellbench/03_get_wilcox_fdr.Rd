mtd = commandArgs(trailingOnly = T)[1]
print(mtd)
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/diff/wilcox/', mtd,'/res/')
dir.create(rdir, showWarnings = F, recursive = T)
imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/sc_10x_5cl.rds'))
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/sc_10x_5cl/genebycell.rds')
raw = raw[,colnames(imp)]
g = intersect(rownames(imp), rownames(raw))
imp = imp[g, ]
raw = raw[g, ]
ct = sub('.*:','',colnames(imp))
allf = sub('_diffgene.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/bulkdiff/'))
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

