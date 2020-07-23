ddir = '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/10xcellline/'
bexpr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/jurkat_hex.rds')
v1 = rowMeans(bexpr[,which(colnames(bexpr)=='239T')])
v2 = rowMeans(bexpr[,which(colnames(bexpr)=='Jurkat')])
bexpr = cbind(v1,v2)
colnames(bexpr) = c('293T','jurkat')
allmtd = list.files(ddir)
allmtd = setdiff(allmtd, 'viper')
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/10x/cellline/hg19/norm_genebycell.rds')
rawg = list()
ct = sub('_.*','',colnames(raw))
for (cl in unique(ct)){
  tmp = raw[,ct==cl]
  rawg[[cl]] = rownames(raw[rowMeans(tmp>1)>=0.1,])
}

for (mtd in allmtd){
  if (length(list.files(paste0(ddir,mtd))) != 0){
    sexpr = readRDS(paste0(ddir,mtd,'/hg19.rds'))
    ct = sub('_.*','',colnames(sexpr))
    res <- lapply(colnames(bexpr), function(cl){
      be <- bexpr[,cl]
      be = be[be>1]
      intgene <- intersect(rawg[[cl]],names(be))
      print(length(intgene))
      a = apply(sexpr[intgene, ,drop=F],2,cor,be[intgene],method='spearman')
    })
    names(res) <- colnames(bexpr)
    saveRDS(res,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/perf/useHighExprGene/10xcellline_cor/',mtd,'.rds'))      
  }
}
