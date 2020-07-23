num.expressed.gene = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/code/qc/num.expressed.gene.rds')
prop.mito.read = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/code/qc/prop.mito.read.rds')

pd1 <- unlist(num.expressed.gene)
pd2 <- unlist(prop.mito.read)
pd <- data.frame(data=c(sub('\\..*','',names(pd1)),sub('\\..*','',names(pd2))),value=c(pd1,pd2),type=c(rep('number.expressed.gene',length(pd1)),rep('proportion.mitochondrial.read',length(pd2))),stringsAsFactors = F)

library(ggplot2)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/code/qc/qc.pdf',height=5.5,width=6.5)
ggplot(pd,aes(x=data,y=value,col=data)) + geom_violin() + geom_jitter(width = 0.1,size=0.1,alpha=0.05) + theme_classic() + facet_wrap(~type,nrow=1,scales = 'free') + coord_flip() + theme(legend.position = 'none') 
dev.off()

