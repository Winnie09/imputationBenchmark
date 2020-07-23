bexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/jurkat_hex.rds')
v1 = rowMeans(bexpr[,which(colnames(bexpr)=='239T')])
v2 = rowMeans(bexpr[,which(colnames(bexpr)=='Jurkat')])
bexpr = cbind(v1,v2)
colnames(bexpr) = c('293T','jurkat')
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/10xcellline/')
for (mtd in allmtd){
      if (length(list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/10xcellline/',mtd))) != 0){
            sexpr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/10xcellline/',mtd,'/hg19.rds'))
            scl = sub('_.*','',colnames(sexpr))
            scl = sub('X293T','293T',scl)
            res <- lapply(colnames(bexpr), function(cl){
                  be <- bexpr[,cl]
                  intgene <- intersect(row.names(sexpr),names(be))
                  apply(sexpr[intgene,which(scl == cl)],2,cor,be[intgene],method='spearman')
            })
            names(res) <- colnames(bexpr)
            saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/10xcellline_cor/',mtd,'.rds'))      
      }
      
}
