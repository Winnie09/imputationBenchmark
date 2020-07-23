trajmethod = commandArgs(trailingOnly = T)[1]
suppressMessages(library('monocle'))
library(TSCAN)
library(ggplot2)
library(gridExtra)
allmtd = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/',trajmethod,'/cds'))
res <- lapply(allmtd,function(mtd){
  print(mtd)
  dir.create(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/plot/plot/',trajmethod, '/cds/'),showWarnings = F,recursive = T)
  pdf(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/plot/plot/',trajmethod, '/cds/',mtd,'.pdf'),width = 10, height = 12)
  allf = list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/',trajmethod,'/cds/',mtd))
  plist = list()
  for (f in allf){
    if (trajmethod=='monocle2'){
      cds = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/trajectory/result/',trajmethod,'/cds/',mtd,'/',f))
      plist[[f]] <- plot_cell_trajectory(cds, cell_size = 0.2) + ggtitle(sub('.rds','',f)) + theme(legend.position = 'none')  
    } else if (trajmethod=='tscan'){
      mcl <- cds[['mcl']]
      row.names(mcl[[1]]) <- paste0(sub(':','_',row.names(mcl[[1]])),1:nrow(mcl[[1]]))
      names(mcl[[3]]) <- paste0(sub(':','_',names(mcl[[3]])),1:nrow(mcl[[1]]))
      plist[[f]] <- plotmclust(mcl) + theme(legend.position = 'none') + ggtitle(sub('.rds','',f))  
    }
  }
  grid.arrange(grobs=plist)
  dev.off()
})
