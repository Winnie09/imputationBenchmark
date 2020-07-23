realfdr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plotdata/realfdr.rds')
reportfdr = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plotdata/reportfdr.rds')

library(ggplot2)
plist = list()
for (mtd in names(realfdr)){
      tmp = realfdr[[mtd]]
      q = sub('.rds','', sub('.*_','', names(tmp)))
      q = factor(q, levels = as.character(sort(unique(as.numeric(q)))))
      pd = data.frame(value = tmp, q = q)
      plist[[mtd]] <- ggplot(data = pd, aes(x=q, y = value, color=q)) + geom_point(size=.5) + theme_classic() + ylab('real fdr') + ggtitle(mtd) + theme(legend.position = 'none')
}

library(gridExtra)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plot/realfdr_by_q_scatter.pdf',width=12,height=8)
grid.arrange(grobs = plist)
dev.off()

plist = list()
for (mtd in names(reportfdr)){
      tmp = reportfdr[[mtd]]
      q = sub('.rds','', sub('.*_','', names(tmp)))
      q = factor(q, levels = as.character(sort(unique(as.numeric(q)))))
      pd = data.frame(value = tmp, q = q)
      plist[[mtd]] <- ggplot(data = pd, aes(x=q, y = value, color=q)) + geom_point(size=.5) + theme_classic() + ylab('claimed fdr') + ggtitle(mtd) + theme(legend.position = 'none')
}

library(gridExtra)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/simu/plot/plot/claimfdr_by_q_scatter.pdf',width=12,height=8)
grid.arrange(grobs = plist)
dev.off()

