## /home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/check_cellline_differential/code/10x.R
numcell1 = as.numeric(commandArgs(trailingOnly = T)[1])
numcell2 = as.numeric(commandArgs(trailingOnly = T)[2])
seed = as.numeric(commandArgs(trailingOnly = T)[3])
mat = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/count/293t_full_genebycell_count.rds')
## subset cells
set.seed(seed)
sid = sample(ncol(mat))
mat = mat[, sid[1: (numcell1 + numcell2)]]
## normailize
libsize = colSums(mat)
mat <- sweep(mat,2,libsize,'/')
mat <- log2(mat + 1)
## remove low expr genes
kid = which(rowMeans(mat > 0) >=0.01)
mat = mat[kid,]
## wilcoxon
pval = sapply(1:nrow(mat), function(i) {
  a = mat[i, (1: numcell1)]
  b = mat[i, ((numcell1+1): (numcell1 + numcell2))]
  names(a) <- names(b) <- NULL
  if (identical(a,b)) {
    1
  } else {
    res = wilcox.test(a, b)
    res$p.value  
  }
})
fdr = p.adjust(pval, method='fdr')
a = c(sum(fdr < 0.05), length(fdr))
saveRDS(a, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/check_cellline_differential/result/fdr/diffNumCell/',numcell1,'_',numcell2,'_',seed,'.rds'))
