library(ggplot2)
library(gridExtra)
pdir1 <- '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/plot/pca/'
pdir2 <- '/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/plot/umap/'
allmtd = readLines('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/resource/impute_method_latent.txt')
af <- list.files(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/raw'))
af <- af[!grepl('mix', af)]

for (f in af){
  print(f)
  noImp <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/raw/', f))
  ct <- sapply(colnames(noImp), function(i) sub('.*:', '', i))
  for (mtd in allmtd){
      print(mtd)
      dir.create(paste0(pdir1, mtd))
      dir.create(paste0(pdir2, mtd))
      
      if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/pd/pc/',mtd,'/', f))){
        if (!grepl('latent', mtd)){
          d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/pd/pc/',mtd,'/', f))
        } else {
          d = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/',f))
          d <- t(d)
        }
        d <- d[names(ct), ]
        u <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/plot/cellbench/pd/umap/',mtd,'/',f))
        u <- u[names(ct), ]
        if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/cellbench/kmeans/', mtd, '/', f))){
          clu <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/cellbench/kmeans/', mtd, '/', f))
          clu <- clu[names(ct)]
          pd <- rbind(data.frame(x = d[,1], y = d[,2], col = as.factor(ct), type = 'celltype'),
                    data.frame(x = d[,1], y = d[,2], col = as.factor(clu), type = 'kmeans'))
          pd2 <- rbind(data.frame(x = u[,1], y = u[,2], col = as.factor(ct), type = 'celltype'),
                    data.frame(x = u[,1], y = u[,2], col = as.factor(clu), type = 'kmeans'))
        } 
        if (file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/cellbench/louvein/', mtd, '/', f))){
          clu <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/clustering/result/cellbench/louvein/', mtd, '/', f))
          clu <- clu[names(ct)]
          pd <- rbind(pd,
                    data.frame(x = d[,1], y = d[,2], col = as.factor(clu), type = 'louvain'))
          pd2 <- rbind(pd2,
                    data.frame(x = u[,1], y = u[,2], col = as.factor(clu), type = 'louvain'))
        }
        p <- ggplot() + geom_point(data=pd, aes(x = x, y = y, col = col), size = 0.05, alpha = 0.8) + 
          guides(color=guide_legend(override.aes = list(size=5, alpha=1))) +
          theme_void() + xlab('PC1') + ylab('PC2') +
          theme(legend.position = 'none') +
          ggtitle(mtd) +
          facet_wrap(~type)
        ggsave(paste0(pdir1, mtd, '/', sub('rds','png',f)), p, width = 10.5,height = 2, dpi = 300)
        
        p <- ggplot() + geom_point(data=pd2, aes(x = x, y = y, col = col), size = 0.05, alpha = 0.8) + 
          guides(color=guide_legend(override.aes = list(size=5, alpha=1))) +
          theme_void() + xlab('UMAP1') + ylab('UMAP2') +
          theme(legend.position = 'none') +
          ggtitle(mtd) +
          facet_wrap(~type)
        ggsave(paste0(pdir2, mtd, '/', sub('rds','png',f)), p, width = 7.5,height = 2.5, dpi = 300)
        rm(d)
      }
  }
} 


