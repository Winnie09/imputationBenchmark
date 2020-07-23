bexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/hm_cellline_combineEncsr.rds')
colnames(bexpr)[which(colnames(bexpr) == 'H1-hESC')] = 'H1'
colnames(bexpr)[which(colnames(bexpr) == 'IMR-90')] = 'IMR90'
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/')
allmtd = setdiff(allmtd,'screcover')
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/GSE81861/GSE81861_Cell_Line_COUNT/norm_genebycell.rds')
ct = sub('_.*','',colnames(raw))
rawg = list()
for (cl in unique(ct)){
  tmp = raw[,ct==cl]
  rawg[[cl]] = rownames(raw[rowMeans(tmp>1)>=0.1,])
}

for (mtd in allmtd){
  sexpr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/',mtd,'/GSE81861_Cell_Line_COUNT.rds'))
  ct = sub('_.*','',colnames(sexpr))
  res <- lapply(colnames(bexpr), function(cl){
    be <- bexpr[,cl]
    be = be[be>1]
    intgene <- intersect(rawg[[cl]],names(be))
    print(length(intgene))
    a = apply(sexpr[intgene, ,drop=F],2,cor,be[intgene],method='spearman')
  })
  names(res) <- colnames(bexpr)
  dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/useHighExprGene/hm_cellline_cor/', showWarnings = F)
  saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/useHighExprGene/hm_cellline_cor/',mtd,'.rds'))
}
