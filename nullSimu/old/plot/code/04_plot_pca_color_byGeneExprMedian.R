library(ggplot2)
library(gridExtra)
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/')
mtd = allmtd[1]
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',mtd))
for (f in allf){
      if (!f == 'cell1k.rds'){
            lib = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/', sub('.rds','',f), '/lib.rds'))      
      } else {
            lib = rep(1,1000)
      }
      for (i in 1:length(allmtd)){
            mtd = allmtd[i]
            print(f)
            lst = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',mtd,'/',f))
            pc = lst$pc
            v = lst$var
            d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/result/procimpute/raw/',f))
            exprMean = apply(d,2,mean)
            if (i == 1) {
                  pd = data.frame(pc1 = pc[,1], pc2 = pc[,2], method = mtd, lib = lib, exprMean = exprMean)
            } else {
                  pd = rbind(pd, data.frame(pc1 = pc[,1], pc2 = pc[,2], method = mtd, lib = lib, exprMean = exprMean))
            }
      }
      pd$method = factor(pd$method, levels = c('raw', setdiff(unique(pd$method),'raw')))
      pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/plot/pca/',sub('.rds','',f),'.pdf'),width=12,height=8)
      if (!f == 'cell1k.rds'){
            print(ggplot(data = pd, aes(x=pc1,y=pc2,color=lib)) + geom_point(size = 0.2) + facet_wrap(~method, scales = 'free') + scale_color_gradient2(low="yellow",mid='green',high="blue",midpoint = mean(pd$lib)) + theme_classic() + theme(legend.position = 'bottom'))      
      } else {
            print(ggplot(data = pd, aes(x=pc1,y=pc2)) + geom_point(size = 0.2) + facet_wrap(~method, scales = 'free') + theme_classic() + theme(legend.position = 'none'))      
      }
      dev.off()      
      
      pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/plot/pca/',sub('.rds','',f),'_colorbyExpr.pdf'),width=14,height=8) 
      print(ggplot(data = pd, aes(x=pc1,y=pc2,color=exprMean)) + geom_point(size = 0.2) + facet_wrap(~method, scales = 'free') + scale_color_gradient2(low="yellow",mid='green',high="blue",midpoint = median(pd$exprMean)) + theme_classic() + theme(legend.position = 'bottom') + guides(fill = guide_colourbar(barwidth = 15)))      
      dev.off()      
}


