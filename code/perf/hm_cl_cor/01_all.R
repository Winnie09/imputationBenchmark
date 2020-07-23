bexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
colnames(bexpr)[which(colnames(bexpr) == 'H1-hESC')] = 'H1'
colnames(bexpr)[which(colnames(bexpr) == 'IMR-90')] = 'IMR90'
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/')
for (mtd in allmtd){
      sexpr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/',mtd,'/GSE81861_Cell_Line_COUNT.rds'))
      ct = sub('_.*','',colnames(sexpr))
      res <- lapply(colnames(bexpr), function(cl){
            be <- bexpr[,cl]
            intgene <- intersect(row.names(sexpr),names(be))
            apply(sexpr[intgene, ct== cl],2,cor,be[intgene],method='spearman')
            
      })
      names(res) <- colnames(bexpr)
      saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/hm_cellline_cor/',mtd,'.rds'))
}
