library(ggplot2)
library(gridExtra)
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/')
allmtd = setdiff(allmtd, 'zinbwave')
allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',allmtd[1]))
for (f in allf){
  if (f != 'cell1k.rds'){
    lib = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/data/processed/', sub('.rds','',f), '/lib.rds'))      
  } else {
    lib = rep(1,1000)
  }
  for (i in 1:length(allmtd)){
    mtd = allmtd[i]
    lst = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/pd/pca/',mtd,'/',f))
    pc = lst$pc
    v = lst$var
    if (i == 1) {
      pd = data.frame(pc1 = pc[,1], pc2 = pc[,2], method = mtd, lib = lib)
    } else {
      pd = rbind(pd, data.frame(pc1 = pc[,1], pc2 = pc[,2], method = mtd, lib = lib))
    }
  }
  pd$method = factor(pd$method, levels = c('raw', setdiff(unique(pd$method),'raw')))
  colnames(pd)[4] = 'LibrarySizeFactor'
  saveRDS(pd,paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/plot/pca/',sub('.rds','',f),'_pd.pdf'))
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/nullSimu/plot/plot/pca/',sub('.rds','',f),'_include_latent.pdf'),width=14,height=5)
  if (!f == 'cell1k.rds'){
    print(ggplot(data = pd, aes(x=pc1,y=pc2,color=LibrarySizeFactor)) + geom_point(size = 0.2) + facet_wrap(~method, scales = 'free',nrow=3) + 
            scale_color_gradient2(low="yellow",mid='green',high="blue",midpoint = mean(pd$LibrarySizeFactor)) + 
            theme_classic() + theme(legend.position = 'bottom')+ xlab('Principal Component 1') + ylab('Principal Component 2')+
            guides(fill = guide_colourbar(barwidth = 8, barheight = 0.5)))
    
  } else {
    print(ggplot(data = pd, aes(x=pc1,y=pc2)) + geom_point(size = 0.2) + facet_wrap(~method, scales = 'free') + theme_classic() + theme(legend.position = 'none'))      
  }
  
  dev.off()      
}

