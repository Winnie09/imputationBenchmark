allf = sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/diff'))
allf = unique(allf)
f = allf[1]
library(ggplot2)
library(gridExtra)
d = list()
for (i in sub('.rds', '', allf)){
      tmp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/result/procdiff/',i,'.rds'))
      d[[i]] = tmp
}

allres <- NULL
for (i in names(d)){
      res <- lapply(d[[i]],function(j) {
            cbind(seq(0,1,0.01),approx(j[,2],j[,1],seq(0,1,0.01),rule=2)$y)
      })
      res <- do.call(rbind,res)
      res <- tapply(res[,2],res[,1],quantile,c(0.25,0.5,0.75))
      allres <- rbind(allres,data.frame(i,as.numeric(names(res)),do.call(rbind,res),stringsAsFactors = F))
}
colnames(allres) <- c('method','x','y1','y2','y3')

pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plot/diff_roc_smooth.pdf'),width=8,height=5)
ggplot(data=allres) +facet_wrap(~method)+ geom_line(data=allres,aes(x=x,y=y2),col='blue') + 
      geom_line(data=allres,aes(x=x,y=y1),col='grey') + geom_line(data=allres,aes(x=x,y=y3),col='grey') + 
      theme_classic() +geom_vline(xintercept = 0.25, colour='red') + xlab('real_fdr') + ylab('sensitivity')
dev.off()

