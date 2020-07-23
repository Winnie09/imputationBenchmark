# ------------
# fluidigm 5cl
# ------------
setwd('/dcl02/hongkai/data/whou')
source('./rna_imputation/resource/function.R')
af <- list.files('./rna_imputation/magnitute/result/hm_cl/')
m <- lapply(af, function(f){
  tmp <- readRDS(paste0('./rna_imputation/magnitute/result/hm_cl/',method,'.rds'))
  data.frame(tmp, method = sub('.rds', '', f), platform = 'fluidigm', stringsAsFactors = FALSE)
})
m_flui <- do.call(rbind, m)

# -------------
# cellbench 5cl
# -------------

af <- list.files('./rna_imputation/magnitute/result/sc_10x_5cl/')
m <- lapply(af, function(f){
  tmp <- readRDS(paste0('./rna_imputation/magnitute/result/sc_10x_5cl/',method,'.rds'))
  data.frame(tmp, method = sub('.rds', '', f), platform = '10x', stringsAsFactors = FALSE)
})
m_umi <- do.call(rbind, m)

# -------
# combine
# -------
pd <- rbind( m_flui, m_umi)
pd[,5] <- factor(pd[,5], levels = c('10x', 'fluidigm'))
pd[,1] <- factor(pd[,1], levels = c("HCC827_H2228", "HCC827_H838", "HCC827_A549", "HCC827_H1975",  "H2228_H838", "H2228_A549",  "H2228_H1975",  "H838_A549","H838_H1975","A549_H1975", "A549_GM12878",  "A549_H1","A549_IMR90","A549_K562", "GM12878_H1","GM12878_IMR90", "GM12878_K562","H1_IMR9", "H1_K562", "IMR90_K562"))
pd[,2] <- as.numeric(pd[,2])
pd[,3] <- as.numeric(pd[,3])



a = tapply(pd$high_LFC_SCC, list(pd$method), mean)
b = tapply(pd$low_LFC_SCC, list(pd$method), mean)
mtdorder_high = names(sort(a))
mtdorder_low = names(sort(b))

library(ggplot2)
library(RColorBrewer)
pd1 = pd
pd1$method <- factor(pd1$method, levels = mtdorder_high)

p_high <- ggplot() + geom_tile(data=pd1, aes(x=celltype, y=method, fill=high_LFC_SCC)) + 
  theme_hm(method_vec=pd1[,'method'])+
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(-0.1,0.0,0.1,0.2,0.3,0.5,0.6,0.7,0.8,1), limits = c(-0.1, 0.9)) + xlab('') + ylab('') + 
  theme(legend.position = 'bottom')+
  labs(fill='Spearman correlation') +
  guides(fill = guide_colourbar(barwidth = 8, barheight = 0.6,title.position = 'top',title.hjust=0.5,vjust=10)) +
  ggtitle('High Log-fold-change')

pd2 = pd
pd2$method <- factor(pd2$method, levels = mtdorder_low)
p_low <- ggplot() + 
  geom_tile(data=pd2, aes(x=celltype, y=method, fill=low_LFC_SCC)) + 
  theme_minimal() +
  theme_hm(method_vec=pd2[,'method'])+
  scale_fill_gradientn(colors=rev(brewer.pal(9,'YlGnBu')),values=c(-0.1,0.0,0.1,0.2,0.3,0.5,0.6,0.7,0.8,1), limits = c(-0.1, 0.9)) + xlab('') + ylab('') + 
  theme(legend.position = 'bottom')+
  labs(fill='Spearman correlation') +
  guides(fill = guide_colourbar(barwidth = 8, barheight = 0.6,title.position = 'top',title.hjust=0.5,vjust=10)) +
  ggtitle('Low Log-fold-change')

library(gridExtra)
pdf('./rna_imputation/magnitute/plot/high_low_lfc_heatmap.pdf', width = 8, height= 5.5)
grid.arrange(p_high, p_low, nrow = 1)
dev.off()

