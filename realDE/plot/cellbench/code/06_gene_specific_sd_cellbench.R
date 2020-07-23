coldf = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/resource/method_latent_color.rds')
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/sc_10x_5cl/genebycell.rds')
ct = sub('.*:','',colnames(raw))

apd <- lapply(unique(ct), function(sct){
  print(sct)
  sd2 = apply(raw[,ct==sct],1,sd)
  names(sd2) = rep('no_imp',length(sd2))
  allmtd = list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/')
  sd_imp <- sapply(allmtd, function(mtd){
    print(mtd)
    imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/sc_10x_5cl.rds'))
    raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/sc_10x_5cl/genebycell.rds')
    imp = imp[,colnames(raw)]
    ct = sub('.*:','',colnames(imp))
    apply(imp[,ct==sct],1,sd)
  })
  
  sd1 = unlist(sd_imp)
  names(sd1) = sub('\\..*','',names(sd1))
  names(sd1) = coldf[match(names(sd1), coldf[,'shortName']),'fullName']  
  pd = data.frame(sd = c(sd1,sd2), method = c(names(sd1), names(sd2)), celltype=sct)
})

saveRDS(apd,'/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/sd_pd.rds')

pd = do.call(rbind, apd)
library(ggplot2)
library(reshape2)
colv = coldf[match(as.character(pd$method), coldf[,'fullName']),'color']
names(colv) = as.character(pd$method)

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/plot/cellbench/plot/sd.pdf',width=10,height=8)
ggplot() + geom_boxplot(data=pd, aes(x=method, y = sd, fill=method),alpha=0.9,outlier.shape = NA) +
  theme_classic() + ylim(c(0,2.5))+ylab('Gene-specific standard deviation')+
  theme(axis.text.x = element_text(angle=90,hjust=1,color='black'), axis.text.y=element_text(color='black'),legend.position = 'none') +
  scale_fill_manual(values=colv) +
  facet_wrap(~celltype)
dev.off()
