clu <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cell_celltype_cluster.rds')
ct <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/cluster/cluster_celltype.rds')
ct <- ct[!is.na(ct[,2]),]

mtd = commandArgs(trailingOnly = T)[1]
print(mtd)
suppressMessages(library(MAST))
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/mast/', mtd,'/res/'), showWarnings = F, recursive = T)
imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/hca/',mtd,'/MantonBM6.rds'))
if (grepl('A.1',colnames(imp)[1])){
  colnames(imp) = sub('\\.','-',colnames(imp))
}
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/hca/MantonBM6/genebycell.rds')
raw = raw[,colnames(imp)]

cid <- expand.grid(ct[,1],ct[,1])
cid <- cid[cid[,1] > cid[,2],]
cid <- cid[ct[match(cid[,1],ct[,1]),2]!=ct[match(cid[,2],ct[,1]),2],]
for (i in 1:nrow(cid)){
  clu1 <- cid[i,1]
  clu2 <- cid[i,2]
  if (!file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/mast/', mtd,'/res/',ct[ct[,1]==clu1,2],'-',clu1,'_',ct[ct[,1]==clu2,2],'-',clu2,'.rds'))){
    count = cbind(raw[,clu[clu[,2]==clu1,1]], raw[,clu[clu[,2]==clu2,1]])
    expr = cbind(imp[,clu[clu[,2]==clu1,1]], imp[,clu[clu[,2]==clu2,1]])
    g = intersect(rownames(expr), rownames(count)) ##
    expr = expr[g,] ##
    count = count[g,] ##
    cdr <- scale(colMeans(count > 0))
    cluster <- rep('clu1',ncol(expr)) ## length #cells
    cluster[(sum(clu[,2]==clu1)+1):ncol(expr)] <- 'clu2' ###
    sca <- FromMatrix(expr,data.frame(wellKey=colnames(expr),cluster=cluster,cngeneson=cdr), data.frame(primerid=row.names(expr),Gene=row.names(expr)), check_sanity = F)
    zlmCond <- zlm(~cluster+cngeneson, sca)
    summaryCond <- summary(zlmCond, doLRT="clusterclu2")
    summaryDt <- as.data.frame(summaryCond$datatable)
    pval <- summaryDt[summaryDt$component == "H" & summaryDt$contrast == "clusterclu2",c(1,4)]
    lfc <- summaryDt[summaryDt$component == "logFC" & summaryDt$contrast == "clusterclu2",c(1,7)]
    combine <- merge(pval,lfc)
    colnames(combine) <- c("Gene","pvalue","Log-foldchange")
    row.names(combine) <- combine[,1]
    res <- combine
    colnames(res)[colnames(res)=='Log-foldchange'] <- "stat"
    res <- res[order(res[,'pvalue'],-abs(res[,'stat'])),]
    saveRDS(res, paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/hca/diff/mast/', mtd,'/res/',ct[ct[,1]==clu1,2],'-',clu1,'_',ct[ct[,1]==clu2,2],'-',clu2,'.rds'))    
  }
}

