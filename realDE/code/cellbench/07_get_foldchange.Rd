mtd = commandArgs(trailingOnly = T)[1]
print(mtd)
rdir = paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/diff/foldchange/', mtd,'/')
dir.create(rdir, showWarnings = F, recursive = T)
imp = readRDS(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/result/procimpute/cellbench/',mtd,'/sc_10x_5cl.rds'))
raw = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/data/processed/cellbench/sc_10x_5cl/genebycell.rds')
raw = raw[,colnames(imp)]
g = intersect(rownames(imp), rownames(raw))
imp = imp[g, ]
raw = raw[g, ]
ct = sub('.*:','',colnames(imp))
allf = sub('_diffgene.rds','',list.files('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/realDE/result/bulkdiff/'))
for (f in allf){
  ct1 = sub('_.*','',f)
  ct2 = sub('.*_','',f)
  fc = rowMeans(imp[,ct==ct1]) - rowMeans(imp[,ct==ct2])
  fc = sort(abs(fc),decreasing = T)
  df = data.frame(Gene=names(fc), fc = fc)
  saveRDS(df, paste0(rdir,f,'.rds'))  
}

