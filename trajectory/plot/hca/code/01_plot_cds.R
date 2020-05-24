trajmethod = commandArgs(trailingOnly = T)[1]
suppressMessages(library('monocle'))
library(TSCAN)
library(ggplot2)
library(gridExtra)
allmtd = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/hca/',trajmethod,'/cds'))
dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/plot/hca/plot/',trajmethod, '/cds/'),showWarnings = F,recursive = T)
pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/plot/hca/plot/',trajmethod, '/cds/all.pdf'),width = 30, height = 30)
plist = list()
for (mtd in allmtd){
  print(mtd)
  f = 'MantonBM6.rds'
  if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/hca/',trajmethod,'/cds/',mtd,'/',f))){
    cds = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/hca/',trajmethod,'/cds/',mtd,'/',f))
    if (trajmethod=='monocle2'){
      
      plist[[mtd]] <- plot_cell_trajectory(cds, cell_size = 0.2) + ggtitle(mtd) + theme(legend.position = 'none')  
    } else if (trajmethod=='tscan'){
      mcl <- cds[['mcl']]
      row.names(mcl[[1]]) <- paste0(sub(':','_',row.names(mcl[[1]])),1:nrow(mcl[[1]]))
      names(mcl[[3]]) <- paste0(sub(':','_',names(mcl[[3]])),1:nrow(mcl[[1]]))
      plist[[mtd]] <- plotmclust(mcl) + theme(legend.position = 'none') + ggtitle(mtd)
    }  
  }
  
}
grid.arrange(grobs=plist)
dev.off()  

