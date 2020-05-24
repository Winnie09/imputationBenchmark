mtd = commandArgs(trailingOnly = T)[1]
print(mtd)
suppressMessages(library(MAST))
library(Matrix)
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullDE/result/10x/diff/mast/', mtd,'/res/')
dir.create(rdir, showWarnings = F, recursive = T)
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',mtd))
q = sapply(allf, function(i) strsplit(sub('.rds','',i), '_')[[1]][4] )
allf = allf[q==0]
a <- sapply(allf, function(f){
  print(f)
  print(f)
  cn1 = as.numeric(strsplit(sub('.rds','',f), '_')[[1]][1])
  cn2 = as.numeric(strsplit(sub('.rds','',f), '_')[[1]][2])
  count = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/', sub('.rds','',f), '/genebycell.rds'))
  expr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',mtd,'/',f))
  colnames(expr) = colnames(count)
  rownames(expr) = rownames(count)
  print(rankMatrix(expr)[1])
  # if (rankMatrix(expr)[1] != ncol(expr)){
  #   set.seed(12345)
  #   u = sample(ncol(expr),rankMatrix(expr)[1])
  #   expr[,u] = expr[,u] + rnorm(nrow(expr))*0.001
  # }
  # print(rankMatrix(expr)[1])
  cdr <- scale(colMeans(count > 0))
  cluster <- rep('clu1',ncol(expr)) ## length = #cells
  cluster[(cn1+1):ncol(expr)] <- 'clu2' ###
  sca <- FromMatrix(expr,data.frame(wellKey=colnames(expr),cluster=cluster,cngeneson=cdr), data.frame(primerid=row.names(expr),Gene=row.names(expr)), check_sanity = FALSE)
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
  saveRDS(res, paste0(rdir,f))  
  return(0)
})
  

