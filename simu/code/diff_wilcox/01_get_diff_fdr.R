library(parallel)
allm <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/')
m = as.character(commandArgs(trailingOnly = T)[1])
print(m)
allf <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m))
res <- lapply(allf, function(f) {
      print(f)
      p = strsplit(f,'_')[[1]][3]
      q = sub('.rds','',strsplit(f,'_')[[1]][4])
      if (p != 0 & q !=0){
            d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procimpute/',m,'/',f))
            pval <- sapply(1:nrow(d), function(i) {
                  if (length(unique(d[i,]))==1){
                        1
                  } else {
                        wilcox.test(d[i,1:(ncol(d)/2)],d[i,(1+(ncol(d)/2)):ncol(d)])$p.value  
                  }
            })
            fdrv <- p.adjust(pval,method='fdr')
            g <- row.names(d)
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
            dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_wilcox/diff/',m), showWarnings = F)
            saveRDS(perf,file=paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_wilcox/diff/',m,'/',f))
      }
})


