# Title     : Run Saver1.0
# Objective : Compare with Saver1.0
# Created by: rui
# Created on: 5/9/18

devtools::install_github("mohuangx/SAVER@v1.0.0")
library('SAVER')
getwd()
packageVersion('SAVER')
sessionInfo()  # see version

# read data (COUNT)
in_file = '~/data/cell_row/mouse_brain.g28k_c1.3m.csv'
out_file = 'mouse_brain.g28k_c1.3m.saver.csv.gz'
df = read.csv(in_file, row.names=1)
df = t(df)  # into gene_row
df[0:3, 0:3]

# single core
# saver4 <- saver(df)

# run in parallel
library(doParallel)
cl <- makeCluster(4, outfile = "")
registerDoParallel(cl)
saver5 <- saver(df)
gc(verbose=TRUE) #RAM usage
stopCluster(cl)


# save result
z <- gzfile(out_file)
write.csv(saver5$estimate, z)