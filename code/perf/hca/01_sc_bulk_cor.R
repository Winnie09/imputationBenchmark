bexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/hca/GSE74246_RNAseq_normalized.rds')
tmpct =  sub('.*\\.','',colnames(bexpr))
bexpr = sapply(unique(tmpct), function(sct) rowMeans(bexpr[,tmpct==sct],na.rm = T))
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/')
allmtd = allmtd[!grepl('latent',allmtd)] 
for (mtd in allmtd){
      sexpr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd,'/MantonBM6.rds'))
      ct = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/ct/MantonBM6.rds')
      ct = ct[match(colnames(sexpr),names(ct))]
      res <- lapply(colnames(bexpr), function(cl){
            be <- bexpr[,cl]
            intgene <- intersect(row.names(sexpr),names(be))
            apply(sexpr[intgene, ct== cl],2,cor,be[intgene],method='spearman')
      })
      names(res) <- colnames(bexpr)
      
      dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/hca/')
      saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/hca/',mtd,'.rds'))
}

