mtd = commandArgs(trailingOnly = T)[1]
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/10x/diff/wilcox/', mtd,'/res/')
dir.create(rdir, showWarnings = F, recursive = T)
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',mtd))
q = sapply(allf, function(i) strsplit(sub('.rds','',i), '_')[[1]][4] )
allf = allf[q==0]
for (f in allf){
  cn1 = as.numeric(strsplit(sub('.rds','',f), '_')[[1]][1])
  cn2 = as.numeric(strsplit(sub('.rds','',f), '_')[[1]][2])
  imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',mtd,'/',f))
  count = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/', sub('.rds','',f), '/genebycell.rds'))
  rownames(imp) = rownames(count)
  colnames(imp) = colnames(count)
  # g = intersect(rownames(expr), rownames(count)) ##
  # expr = expr[g,] ##
  # count = count[g,] ##
  pval <- sapply(1:nrow(imp), function(i) {
        v1 = imp[i, 1:cn1]
        v2 = imp[i, ((cn1+1):ncol(imp))]
        if (length(unique(c(v1,v2)))==1){
              1
        } else {
              wilcox.test(v1, v2)$p.value  
        }
  })
  fdrv <- p.adjust(pval,method='fdr')
  df = data.frame(Gene=rownames(imp), fdr = fdrv)
  saveRDS(df, paste0(rdir,f))  
}

