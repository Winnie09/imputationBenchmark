rm(list=ls())
#fixed parameters for GiniClust2
minCellNum           = 3 # filtering, remove genes expressed in fewer than minCellNum cells
minGeneNum           = 200 # filtering, remove cells expressed in fewer than minGeneNum genes
expressed_cutoff     = 0 # filtering, for raw counts
gini.bi              = 0 # fitting, default is 0, for qPCR data, set as 1. 
log2.expr.cutoffl    = 0 # cutoff for range of gene expression   
log2.expr.cutoffh    = 20 # cutoff for range of gene expression 
Gini.pvalue_cutoff   = 0.05 # fitting, Pvalue, control how many Gini genes chosen
Norm.Gini.cutoff     = 1 # fitting, NormGini, control how many Gini genes chosen, 1 means not used.
span                 = 0.9 # parameter for LOESS fitting
outlier_remove       = 0.75 # parameter for LOESS fitting
#GeneList             = 1  # parameter for clustering, 1 means using pvalue, 0 means using HighNormGini
#Gamma                = 0.9 # parameter for clustering
diff.cutoff          = 1# MAST analysis, filter genes that don't have high log2_foldchange to reduce gene num
lr.p_value_cutoff    = 1e-5 # MAST analysis, pvalue cutoff to identify differentially expressed genes
CountsForNormalized  = 10000 # if normalizing- by default not used

#paths
args <- commandArgs()
baseName <- args[6]
Rfundir <- paste(baseName, "GiniClust2/", sep="/")
workdir <- Rfundir

exprimentID          = "experiment_name"                                	# experiment or data set ID

#setwd(workdir)
dir.create(file.path(workdir, "results"), showWarnings = FALSE) #folder to save results
dir.create(file.path(workdir, "figures"), showWarnings = FALSE) #folder to save figures
#load packages and functions
source(paste(Rfundir,"GiniClust2_packages.R",sep=""))
source(paste(Rfundir,"GiniClust2_functions.R",sep=""))


# infile1 <- paste(getwd(), baseName, "raw.txt", sep="/")
infile1 <- paste(getwd(), "raw_temporal.txt", sep="/") 
data <- read.table(file=infile1, sep = '\t',row.names=1,header=T)


source(paste(Rfundir,"GiniClust2_preprocess.R",sep=""))
source(paste(Rfundir,"GiniClust2_filtering_RawCounts.R",sep=""))

#Gini-based clustering steps
source(paste(Rfundir,"GiniClust2_fitting.R",sep=""))
#load(paste("GiniClust2/results/", exprimentID,"_StatisticsTable_afterLOESS.RData",sep=""))  #"Genelist.top_pvalue" "ExprM.Stat2" 
#load(paste("GiniClust2/results/", exprimentID,"_ExprM.filter.RData", sep=""))  #"ExprM.RawCounts.filter"  "ExprM.normCounts.filter"           
#load(paste("GiniClust2/results/", exprimentID,"_ExprM.RData", sep=""))  #"ExprM.RawCounts"  "ExprM.normCounts"           

#dim(ExprM.RawCounts.filter)

GeneList.final = Genelist.top_pvalue
Genes <- rownames(data)
id<- match(GeneList.final, Genes)

out1 <- paste(baseName,  "Gini_ID.txt", sep="/")
write.table(file= out1, id, sep="\t", row.names = FALSE, col.names = FALSE)
