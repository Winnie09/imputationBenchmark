#Preprocess data

#raw data
ExprM.RawCounts <- data

if(length(colnames(ExprM.RawCounts))==0){
  colnames(ExprM.RawCounts)<-paste("C",1:dim(ExprM.RawCounts)[2],sep="")
}

if(length(rownames(ExprM.RawCounts))==0){
  print("Please label rows with gene names.")
}
#normalized data
colsum <- apply(ExprM.RawCounts,2,sum)
ExprM.normCounts <- t(t(ExprM.RawCounts)*CountsForNormalized/colsum)

write.table(ExprM.RawCounts, file=paste("GiniClust2/results/", exprimentID, "_rawCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(ExprM.normCounts, file=paste("GiniClust2/results/", exprimentID, "_normCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
save(ExprM.RawCounts, ExprM.normCounts, file=paste("GiniClust2/results/",exprimentID, "_ExprM.RData",sep=""))
