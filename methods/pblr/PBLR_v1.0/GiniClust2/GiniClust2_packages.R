#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("hydroGOF")) {
  install.packages("hydroGOF", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("glmnet")) {
  install.packages("glmnet", dependencies = TRUE, repos="http://cran.r-project.org")
}
#if (!require("doMC")) {
#  install.packages("doMC", dependencies = TRUE, repos="http://cran.r-project.org")
#}
if (!require("StatMatch")) {
  install.packages("StatMatch", dependencies = TRUE, repos="http://cran.r-project.org")
}
# if (!require("Rtsne")) {
#   install.packages("Rtsne", dependencies = TRUE, repos="http://cran.r-project.org")
# }
if (!require("fpc")) {
  install.packages("fpc", dependencies = TRUE, repos="http://cran.r-project.org")
}
# if (!require("GA")) {
#   install.packages("GA", dependencies = TRUE, repos="http://cran.r-project.org")
# }
if (!require("MASS")) {
  install.packages("MASS", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("session")) {
  install.packages("session", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("Matrix")) {
  install.packages("Matrix", dependencies = TRUE, repos="http://cran.r-project.org")
}
# if (!require("vegan")) {
#   install.packages("vegan", dependencies = TRUE, repos="http://cran.r-project.org")
# }
if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("reshape")) {
  install.packages("reshape", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("abind")) {
  install.packages("abind", dependencies = TRUE, repos="http://cran.r-project.org")
}
if (!require("drc")) {
  install.packages("drc", dependencies = TRUE, repos="http://cran.r-project.org")
}
# if (!require("cluster")) {
#   install.packages("cluster", dependencies = TRUE, repos="http://cran.r-project.org")
# }
# if (!require("R.utils")) {
#   install.packages("R.utils", dependencies = TRUE, repos="http://cran.r-project.org")
# }
# if (!require("MAST")) {
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("MAST")
# }
# if (!require("svd")) {
#   install.packages("svd", dependencies = TRUE, repos="http://cran.r-project.org")
# }
# if (!require("GraphAT")) {
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("GraphAT")
# }
# if (!require("varhandle")) {
#   install.packages("varhandle", dependencies = TRUE, repos="http://cran.r-project.org")
# }
# if (!require("gridExtra")) {
#   install.packages("gridExtra", dependencies = TRUE, repos="http://cran.r-project.org")
# }
# if (!require("VennDiagram")) {
#   install.packages("VennDiagram", dependencies = TRUE, repos="http://cran.r-project.org")
# }
# if (!require("scatterplot3d")) {
#   install.packages("scatterplot3d", dependencies = TRUE, repos="http://cran.r-project.org")
# }
# if (!require("RColorBrewer")) {
#   install.packages("RColorBrewer", dependencies = TRUE, repos="http://cran.r-project.org")
# }
if(!require("ROCR")) { install.packages("ROCR"); require("ROCR") }

#Load libraries
library(ggplot2)
library(hydroGOF)
library(glmnet)
#library(doMC)
library(StatMatch)
# library(Rtsne)
library(fpc)
#library(GA)
library(MASS)
library(session)
library(Matrix)
# library(vegan)
library(data.table)
library(plyr)
library(reshape)
library(abind)
library(drc)
# library(cluster)
# library(R.utils)
# library(MAST)
#library(svd)
# library(GraphAT)
# library(varhandle)
# library(gridExtra)
# library(VennDiagram)
# library(scatterplot3d)
# library(RColorBrewer)

#extra RaceID packages
#require(tsne)
#require(pheatmap)
#require(cluster)
#require(mclust)
#require(flexmix)
#require(lattice)
#require(amap)
#require(RColorBrewer)
#require(locfit)
#require(methods)
