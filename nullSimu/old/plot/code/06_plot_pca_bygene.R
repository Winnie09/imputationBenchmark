library(ggplot2)
library(gridExtra)
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/')
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',allmtd[1]))
for (f in allf){
      pd = NULL
      for (i in 1:length(allmtd)){
            mtd = allmtd[i]
            print(f)
            lst = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca_bygene/',mtd,'/',f))
            pc = lst$pc
            pd = rbind(pd,data.frame(method = mtd,pc1 = pc[,1], pc2 = pc[,2], geneExprMean = lst$geneExprMean, geneExprVar = lst$geneExprVar, geneLen = lst$geneLen))
      }
      pd$method = factor(pd$method, levels = c('raw', setdiff(unique(pd$method),'raw')))
      dir.create('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/plot/pca_bygene/', showWarnings = F)
      pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/plot/pca_bygene/',sub('.rds','',f),'.pdf'),width=12,height=8)
      print(ggplot(data = pd, aes(x=pc1,y=pc2,color=geneExprVar)) + geom_point(size = 0.2) + facet_wrap(~method, scales = 'free') + scale_color_gradientn(colors = rainbow(5)) + theme_classic() + theme(legend.position = 'bottom') + labs(colour = "gene var") )
      print(ggplot(data = pd, aes(x=pc1,y=pc2,color=geneExprMean)) + geom_point(size = 0.2) + facet_wrap(~method, scales = 'free') + scale_color_gradientn(colors = rainbow(5)) + theme_classic() + theme(legend.position = 'bottom') + labs(colour = "gene expr") )
      dev.off()     
}

