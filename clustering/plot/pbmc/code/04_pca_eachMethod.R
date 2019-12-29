clumethod = as.character(commandArgs(trailingOnly = T)[[1]])
library(ggplot2)
library(umap)
library(gridExtra)
#########
allmtd = sub('.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/pd/pc/'))
for (mtd in allmtd){
  umap = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/pd/umap/',mtd,'.rds'))
  pd1 = data.frame(pc1 = umap[,1],pc2=umap[,2])
  ct = factor(sub(':.*','', rownames(umap)),levels=c("cd34","cd4_t_helper","b_cells","cytotoxic_t","cd56_nk","naive_t","memory_t","regulatory_t","cd14_monocytes","naive_cytotoxic"))
  p1 <- ggplot(data=pd1, aes(x = pc1, y = pc2, color = ct)) + geom_point(size = .1, alpha=0.1) + theme_classic() + 
    theme(legend.position = 'right', legend.title=element_blank(), legend.text = element_text(size=8), legend.spacing.x = unit(0.01, 'cm'),legend.spacing.y = unit(0.01, 'cm')) + 
    ggtitle(mtd) + xlab('PC1') + ylab('PC2') +
    guides(color = guide_legend(nrow = 4, byrow = F,override.aes = list(size=1,alpha=1)))
  
  clu = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/pbmc/',clumethod,'/',mtd,'/sorted.rds'))
  cluster = factor(clu[rownames(umap)])
  p2 <- ggplot(data=pd1, aes(x = pc1, y = pc2, color = cluster)) + geom_point(size = .1, alpha=0.1) + theme_classic() +
    theme(legend.position = 'right', legend.title=element_text(color='black'), legend.text = element_text(size=8)) + 
    ggtitle(mtd)  + xlab('PC1') + ylab('PC2') +
    guides(color = guide_legend(nrow  = 4, byrow = F,override.aes = list(size=1,alpha=1)))
  
  dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/plot/',clumethod,'/pca/'), showWarnings = F, recursive = T)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/pbmc/plot/',clumethod,'/pca/',mtd,'.pdf'),width=3.5,height=6)
  grid.arrange(p1,p2,layout_matrix=matrix(c(1,2),2))
  dev.off()
}

