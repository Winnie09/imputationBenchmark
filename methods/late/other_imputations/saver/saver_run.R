# terse #

devtools::install_github("mohuangx/SAVER@v0.4.0")
library('SAVER')
packageVersion('SAVER')
sessionInfo()

# setwd('/mnt/lfs2/rui/scImpute_results/1801_PBMC_g5561/msk98/saver')
getwd()

# read data (COUNT)
df = read.csv('10x_human_pbmc_68k.nz40.msk90.csv.gz', row.names=1)
df[0:3, 0:3]

# 10x_hdf5
library(cellrangerRkit)
fname = '~/imputation/data/10x_mouse_brain_1.3M/1.3M/mouse_brain.10kg.h5'
genome = 'mm10'
GeneBCMatrix = get_matrix_from_h5(fname ,genome)
df = exprs(GeneBCMatrix)  # gene_row

# run in parallel
library(doParallel)
registerDoParallel(cores = 16)  # v0.4
saver4 <- saver(df)  # gene_row

# save result
write.csv(saver4$estimate, file='mouse_brain.10kg.saver.csv')









# devtools::install_github("mohuangx/SAVER")

# setwd('/mnt/lfs2/rui/scImpute_results/1801_PBMC_g5561/msk98/saver')

# df = t(df) # if cell_row

# predictions for top 5 highly expressed genes
# saver1 <- saver(df, npred = 5)

# one CPU
# saver3 <- saver(df)

# save result (Gene_row) (Can mess up _ in cell_id)
# write.csv.gz(saver4$estimate, file='10x_human_pbmc_68k.G9987.csv.gz')

# save obj for future use
# save(saver4, file='saver4.rda')