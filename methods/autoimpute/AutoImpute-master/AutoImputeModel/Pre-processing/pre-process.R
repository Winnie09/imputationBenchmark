# Path definition
rm(list=ls())
source_dir = getwd()
setwd(source_dir)
DATA_DIR <- "../data/raw/"
OUTPUT_DIR <- "../data/"

# Input dataset name
args <- commandArgs(TRUE)
Input_file <- paste(DATA_DIR, args[1], sep='')
print(Input_file)
DATASET <- args[2]
Output_file <- paste(OUTPUT_DIR,DATASET,".csv",sep='')

# Import necessary libraries and functions required for preporcessing
suppressMessages(source("functions.R"))
suppressMessages(source("libraries.R"))

# Read raw data
raw <- read.table(file=Input_file, header = TRUE, sep=",", row.names = 1)

# Pre-process
l <- normalize_by_umi_2(t(raw))
processed <- matrix.subset(l,1000)

# Write the processed data to the output file
write.table(as.matrix(processed), file=paste(Output_file,sep=''), quote=FALSE, sep=",", row.names = TRUE, col.names = TRUE)
