DEMethod = commandArgs(trailingOnly = T)[1]
library(reshape2)
library(ggplot2)
# allf = sub('.rds','',list.files('/Volumes/My Passport/wenpin/rna_imputation/simu/result/procdiff/'))
allf = sub('.rds','',list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/')))
allf = setdiff(allf,'zinbwave')
apx <- sapply(allf, function(f){ ###############3
      print(f)
      # res = readRDS(paste0('/Volumes/My Passport/wenpin/rna_imputation/simu/result/procdiff/',f,'.rds'))
      res = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff_',DEMethod,'/procdiff/',f,'.rds'))
      tmp <- sapply(1:length(res), function(i) { ###########
            i <- res[[i]]
            y = i[,"Real_FDR"]
            x = i[,'Reported_FDR']
            if (length(unique(x))==1) {
                  NA
            } else {
                  approx(x,y,xout=seq(0.0,0.25,length.out = 1000))$y            
            }
      })
      if (is.list(tmp)){
            tmp = do.call(cbind,tmp)
      }
      rowMeans(tmp)
})

auc <- sapply(colnames(apx), function(i){
      a = apx[, i]
      b = seq(0.0,0.25,length.out = 1000)
      sum((a[-length(a)]+a[-1]) * diff(b)/2)
})
colnames(apx) = paste0(colnames(apx),'(',round(auc[match(colnames(apx), names(auc))],2),')')
#colnames(apx) = colnames(apx)[order(auc)]
pd = melt(apx)
pd$Var1 = rep(seq(0.0,0.25,length.out = 1000), length(unique(pd$Var2)))
colnames(pd) = c('x','method','y')
pd$method = factor(pd$method, levels = colnames(apx)[order(auc)])
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/',DEMethod,'/plot/fdr_alpha.pdf'),height=7,width=8)
ggplot(data=pd,aes(x=x, y=y, color=method)) + geom_line() +
      theme_classic() + theme(legend.position = 'top', legend.title=element_blank()) + 
      ylab('fdr')+xlab(expression(alpha)) +
      geom_abline(slope=1, intercept = 0, size=1.5,linetype=2, alpha=.8)
dev.off()


