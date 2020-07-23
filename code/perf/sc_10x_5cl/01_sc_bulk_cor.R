bexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count_average_replicates.rds')

allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/')
for (mtd in allmtd){
      sexpr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/sc_10x_5cl.rds'))
      ct = sub('.*:','',colnames(sexpr))
      res <- lapply(colnames(bexpr), function(cl){
            be <- bexpr[,cl]
            intgene <- intersect(row.names(sexpr),names(be))
            apply(sexpr[intgene, ct== cl],2,cor,be[intgene],method='spearman')
      })
      names(res) <- colnames(bexpr)
      saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/sc_10x_5cl/',mtd,'.rds'))
}

