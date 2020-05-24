trajmtd = commandArgs(trailingOnly = T)[1]
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
source('./resource/function.R')
colordf = readRDS('./resource/method_latent_color.rds')

### compare plotdata
cor1 = readRDS('./result/perf/assess/trajectory_cor_monocle2.rds')
cor2 = readRDS('./result/perf/assess/trajectory_cor_tscan.rds')
ov1 = readRDS('./result/perf/assess/trajectory_ov_monocle2.rds')
ov2 = readRDS('./result/perf/assess/trajectory_ov_tscan.rds')

library(reshape2)
library(ggplot2)
library(ggrepel)
library(gridExtra)
int <- intersect(names(cor1),names(cor2))
pd1 <- data.frame(monocle2=cor1[int],tscan=cor2[int],mtd=int)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/plot/cellbench/plot/monocle2_tscan_compare.pdf',width=8,height=4)
p1 <- ggplot(pd1,aes(x=monocle2,y=tscan,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + 
  theme_minimal() + theme(legend.position = 'none') + ggtitle('correlation')

int <- intersect(names(ov1), names(ov2))
pd2 <- data.frame(monocle2 = ov1[int], tscan = cor2[int], mtd=int)
p2 <- ggplot(pd2, aes(x=monocle2,y=tscan,label=mtd,color=mtd)) + geom_point() + geom_text_repel() + 
  theme_minimal() + theme(legend.position = 'none') + ggtitle('overlap')
grid.arrange(p1,p2, nrow=1)
dev.off()
