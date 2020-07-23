suppressMessages(library(MAST))
library(reshape2)
## Use MAST to do DE (count, expr)
m = as.character(commandArgs(trailingOnly = T)[1])
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_mast/diff/',m), showWarnings = F)
print(m)
allf <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m))
tmp <- sapply(allf[11:length(allf)], function(f){
  print(f)
  p = strsplit(f,'_')[[1]][3]
  q = sub('.rds','',strsplit(f,'_')[[1]][4])
  if (p != 0 & q !=0){
    # tryCatch({
      expr = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m,'/',f))
      count = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/', sub('.rds','',f),'/genebycell.rds'))
      g = intersect(rownames(expr), rownames(count)) ##
      expr = expr[g,] ##
      count = count[g,] ##
      cdr <- scale(colMeans(count > 0))
      cluster <- rep('clu1',ncol(expr)) ## length #cells
      as.numeric(strsplit(sub('.rds','',f), '_')[[1]][1])+1
      cluster[(as.numeric(strsplit(sub('.rds','',f), '_')[[1]][1])+1):ncol(expr)] <- 'clu2' ###
      sca <- FromMatrix(expr,data.frame(wellKey=colnames(expr),cluster=cluster,cngeneson=cdr), data.frame(primerid=row.names(expr),Gene=row.names(expr)),check_sanity = FALSE)
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
      
      ## get fdr, sensitivity
      fdrv <- p.adjust(res[,'pvalue'],method='fdr')
      g <- row.names(res)
      names(fdrv) <- g
      fdrv <- sort(fdrv)
      dg = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/data/processed/diffNumCell/',sub('.rds','',f),'/diffgn.rds'))
      
      diff <- rep('nodiff',length(g))
      names(diff) <- g
      diff[g %in% dg] <- 'diff'
      sumdifffval <- sum(diff=='diff')
      
      diff <- diff[names(fdrv)]
      perf <- t(sapply(1:length(diff), function(i) {
        num <- sum(diff[1:i]=='diff')
        TN <- sum(diff[(i+1):length(diff)]=='nodiff')
        c(num/sumdifffval,(i-num)/i,fdrv[i],TN/sum(diff=='nodiff'))
      }))
      if (nrow(perf) > 1) {
        for (i in (nrow(perf)):2) {
          if (perf[i-1,2] > perf[i,2]) perf[i-1,2] = perf[i,2]
        }
      }
      colnames(perf) <- c('Sensitivity','Real_FDR','Reported_FDR','Specificity')
      perf <- rbind(c(0,0,0,0),perf)
      print(dim(perf))
      saveRDS(perf,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_mast/diff/',m,'/',f))
    # },error=function(e) {},warings=function(w) {})
  }
  return(0)
})







































:wq
