#########################################################
#Filter data
#########################################################
#parameters: expressed_cutoff #value lower than this value could be just noise.

load(paste("GiniClust2/results/", exprimentID, "_ExprM.RData",sep=""))
ExpressedinCell_per_gene=apply(ExprM.RawCounts,1,function(x) length(x[x > expressed_cutoff ]))
nonMir = grep("MIR|Mir", rownames(ExprM.RawCounts), invert = T)  # because Mir gene is usually not accurate 
Genelist = intersect(rownames(ExprM.RawCounts)[nonMir],rownames(ExprM.RawCounts)[ExpressedinCell_per_gene >= minCellNum])
ExpressedGene_per_cell=apply(ExprM.RawCounts[Genelist,],2,function(x) length(x[x>0]))
length(Genelist)
ExprM.RawCounts.filter = ExprM.RawCounts[Genelist,ExpressedGene_per_cell >= minGeneNum]
ExprM.normCounts.filter  = ExprM.normCounts[rownames(ExprM.RawCounts.filter), colnames(ExprM.RawCounts.filter)]
dim(ExprM.RawCounts.filter) 
dim(ExprM.normCounts.filter) 
write.table(ExprM.RawCounts.filter, file=paste("GiniClust2/results/", exprimentID, "_gene.expression.matrix.normCounts.filtered.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE)
write.table(ExprM.normCounts.filter,  file=paste("GiniClust2/results/", exprimentID, "_gene.expression.matrix.RawCounts.filtered.csv", sep=""),  sep=",", row.names = TRUE, col.names = TRUE)
save(ExprM.RawCounts.filter, ExprM.normCounts.filter, file=paste("GiniClust2/results/", exprimentID, "_ExprM.filter.RData", sep=""))