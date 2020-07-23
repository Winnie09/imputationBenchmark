## diff p val
# method = commandArgs(trailingOnly = T)[1]
method = 'saver'
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/function.R')
bulk = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/bulkrna/expr/hm_cellline_encsr.rds')
colnames(bulk)[which(colnames(bulk)=='H1-hESC')] = 'H1'
colnames(bulk)[which(colnames(bulk)=='IMR-90')] = 'IMR90'
bct = sub('.*_','',colnames(bulk))


raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/raw/GSE81861_Cell_Line_COUNT.rds')
if (method == 'raw'){
  sexpr = raw
} else {
  sexpr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/GSE81861/',method,'/GSE81861_Cell_Line_COUNT.rds'))  
  colnames(sexpr) = colnames(raw)
}
ct = sub('_.*','',colnames(sexpr))
intgene = intersect(rownames(bulk),rownames(sexpr))
bulk = bulk[intgene,]
sexpr = sexpr[intgene,]

## scRNA-seq wilcoxon test:
sct = ct[1]
sct2 = ct[78]
pval <- apply(sexpr,1,function(i) {
  if (sd(i[ct==sct])==0 & sd(i[ct==sct2])==0) {
    1
  } else {
    t.test(i[ct==sct],i[ct==sct2])$p.value  
  }
})
#sfc <- rowMeans(sexpr[,ct==sct]) - rowMeans(sexpr[,ct==sct2])
pval[is.na(pval)] = 1
sfdr = p.adjust(pval,method='fdr')

## bulk limma:
library(limma)
fit <- lmFit(bulk[,bct %in% c(sct,sct2)],design=cbind(1,as.numeric(bct[bct %in% c(sct,sct2)]==sct)))
fit <- eBayes(fit)
bfdr <- topTable(fit,n=nrow(bulk),coef=2)[names(sfdr),]$adj.P.Val


## 
library(pROC)
keepid <- which(rowMeans(sexpr >= 1) >= 0.1)
roc_obj <- roc((bfdr < 0.05)[keepid], sfdr[keepid])
auc(roc_obj)
