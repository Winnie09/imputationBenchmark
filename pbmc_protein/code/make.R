library(ROCR)
am <- list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/')
am <- am[!grepl('latent',am)]
ex <- sapply(am,function(i) file.exists(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/',i,'/sorted.rds')))
am <- am[ex]
perf <- sapply(am,function(impmet) {
  print(impmet)
  d <- readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/pbmc/',impmet,'/sorted.rds'))
  fullct <- sub(':.*','',colnames(d))
  sapply(c('CD19','CD14','CD34','CD3D','CD4','CD8A'),function(gene) {
    if (gene=='CD19') {
      ct <- c('b_cells')
    } else if (gene=='CD14') {
      ct <- 'cd14_monocytes'
    } else if (gene=='CD34') {
      ct <- 'cd34'
    } else if (gene=='NCAM1') {
      ct <- 'cd56_nk'
    } else if (gene=='CD3D') {
      ct <- c('cd4_t_helper','cytotoxic_t','memory_t','naive_cytotoxic','naive_t','regulatory_t')
    } else if (gene=='CD4') {
      ct <- c('cd4_t_helper','memory_t','naive_t','regulatory_t')
    } else if (gene=='CD8A') {
      ct <- c('cytotoxic_t','naive_cytotoxic')
    }
    
    truect <- as.numeric(fullct %in% ct)
    geneexpr <- d[gene,]  
    pred <- prediction(geneexpr, truect)
    perf <- performance(pred,measure='auc')@y.values[[1]]
  })
})

saveRDS(perf,file='/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/pbmc_protein/res/perf.rds')
