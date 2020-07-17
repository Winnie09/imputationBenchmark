list.of.packages <- c("R.matlab", "Matrix", "ggplot2", "Rtsne", "svd", "dplyr", "plyr",
                      "data.table", "mclust", "flexclust", "reshape2",
                      "irlba","dynamicTreeCut", "RColorBrewer","GenSA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# ----------------------------
# load relevant libraries
# ----------------------------
library(R.matlab)
library(Matrix)
library(ggplot2)
library(Rtsne)
library(svd)
library(plyr)
library(dplyr)
library(data.table)
library(mclust)
library(flexclust)
library(reshape2)
library(irlba)
library(dynamicTreeCut)
library(RColorBrewer)
library(GenSA)
library(gplots)
# -------------------------------------
