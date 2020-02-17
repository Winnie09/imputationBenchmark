coldf = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/resource/method_latent_color.rds')
allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/')
sd_imp <- sapply(allmtd, function(mtd){
  imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/sc_10x_5cl.rds'))
  raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/sc_10x_5cl/genebycell.rds')
  imp = imp[,colnames(raw)]
  ct = sub('.*:','',colnames(imp))
  apply(imp[,ct=='A549'],1,sd)
})

sd1 = unlist(sd_imp)
names(sd1) = sub('\\..*','',names(sd1))
names(sd1) = coldf[match(names(sd1), coldf[,'shortName']),'fullName']
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/sc_10x_5cl/genebycell.rds')
ct = sub('.*:','',colnames(raw))
sd2 = apply(raw[,ct=='A549'],1,sd)
names(sd2) = rep('no_imp',length(sd2))

library(ggplot2)
library(reshape2)
pd = data.frame(sd = c(sd1,sd2), method = c(names(sd1), names(sd2)))
colv = coldf[match(as.character(pd$method), coldf[,'fullName']),'color']
names(colv) = as.character(pd$method)

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/A549_sd.pdf',width=5,height=5)
ggplot() + geom_boxplot(data=pd, aes(x=method, y = sd, fill=method),outlier.shape = NA) +
  theme_classic() + ylim(c(0,2.5))+
  theme(axis.text.x = element_text(angle=90,hjust=1), legend.position = 'none') +
  scale_fill_manual(values=colv) 
dev.off()
