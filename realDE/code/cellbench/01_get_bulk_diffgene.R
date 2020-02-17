dataset = commandArgs(trailingOnly = T)[[1]]
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/',dataset,'/bulkdiff/'), showWarnings=F, recursive=T)
if (dataset == 'cellbench'){
  cnt = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/cellbench/GSE86337_processed_count.rds')
  ct = sub('_.*','', colnames(cnt))  
} else if (dataset == 'hca'){
  cnt <- read.table("/home-4/zji4@jhu.edu/scratch/scdata/HCA/geobulk/GSE74246/data/GSE74246_RNAseq_All_Counts.txt",header=T)
  row.names(cnt) <- cnt[,1]
  cnt <- as.matrix(cnt[,-1])
  cnt <- cnt[,grep("^X",colnames(cnt))]
  ct <- sub('.*\\.','',colnames(cnt))
  for (i in unique(ct)){
    colnames(cnt)[ct == i] <- paste0(i,'_',1:sum(ct==i))
  }
}

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
      saveRDS(gs,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/',dataset,'/bulkdiff/',i,'_',j,'_diffgene.rds'))
  }
}

