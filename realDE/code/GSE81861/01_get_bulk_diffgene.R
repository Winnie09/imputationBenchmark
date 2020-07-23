dataset = 'GSE81861'
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/',dataset,'/bulkdiff/'), showWarnings=F, recursive=T)
cnt = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/bulk_of_GSE81861_with_replicates_TPM.rds')
ct = sub('_.*','',colnames(cnt))
for (n1 in 1:(length(unique(ct))-1)){
  i = unique(ct)[n1]
  for (n2 in ((n1+1):length(unique(ct)))){
    j = unique(ct)[n2]
    print(paste0(i,'_',j))
    expr = cnt[, ct %in% c(i,j)]
    expr = expr[rowSums(expr>0)>0,]
    sct <- sub('_.*','',colnames(expr))
    library(limma)
    des <- cbind(1,ifelse(sct==i,1,0))
    fit <- eBayes(lmFit(voom(expr,des),design=des))
    res <- topTable(fit,coef=2,number=nrow(expr))
    res <- res[res[,'adj.P.Val']<0.05,]
    gs <- rownames(res)
    str(gs)
    saveRDS(gs,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/',dataset,'/bulkdiff/',i,'_',j,'_diffgene.rds'))
  }
}
