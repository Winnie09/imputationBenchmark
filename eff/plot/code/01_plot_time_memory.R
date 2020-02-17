setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
colordf = readRDS('./resource/method_latent_color.rds')
tb = readRDS('./eff/result/time.rds')
t = matrix(as.numeric(tb[,-1]),nrow=nrow(tb))
dimnames(t) = list(tb[,1],colnames(tb)[2:5])
library(reshape2)
library(ggplot2)
colnames(t) <- c(1000,5000,50000,100000)
rownames(t) <- colordf[match(rownames(t), colordf$shortName),'fullName']
sca =  readRDS('./result/perf/assess/scalability.rds')
names(sca) = colordf[match(names(sca), colordf$shortName),'fullName']
rownames(t) = paste0(rownames(t),'(',round(sca[rownames(t)],3),')')
timedata <- melt(t)
colnames(timedata) <- c('method','cellnum','time')
timedata$cellnum <- log10(timedata$cellnum)
timedata$col  = colordf[match(as.character(sub('\\(.*','',timedata$method)), colordf$fullName),'color']
timedata$method = factor(as.character(timedata$method),levels= rownames(t)[match(names(sort(sca,decreasing = T)), sub('\\(.*','',rownames(t)))])
colv = timedata$col
names(colv) = as.character(timedata$method)

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/plot/plot/time.pdf',height=5.5,width=5.3)
ggplot()+geom_line(data=timedata,aes(x=cellnum,y=time,color=method),linetype = "dotted") + 
  geom_point(data=timedata,aes(x=cellnum,y=time,color=method))+
  scale_color_manual(values=colv) +
  theme_minimal() + xlab(expression(paste('Number of cells (', "log"[10],'-scaled)'))) + ylab('Time (min)')
dev.off()

tb = readRDS('./eff/result/memory.rds')
tb = sub('K','',tb)
t = matrix(as.numeric(tb[,-1]),nrow=nrow(tb))
dimnames(t) = list(tb[,1],c(1000,5000,50000,100000))
rownames(t) <- colordf[match(rownames(t), colordf$shortName),'fullName']
t = t/(1024^2)
mtdorder = names(sort(t[,2],decreasing=T))
memdata <- melt(t)
colnames(memdata) <- c('method','cellnum','memory')
memdata$cellnum <- log10(memdata$cellnum)
memdata$method = factor(as.character(memdata$method),levels=mtdorder)
memdata$col  = colordf[match(memdata$method, colordf$fullName),'color']
colv2 = memdata$col
names(colv2) = as.character(memdata$method)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/eff/plot/plot/memory.pdf',height=5.5,width=5)
ggplot()+geom_line(data=memdata,aes(x=cellnum,y=memory,color=method),linetype = "dotted") + 
  geom_point(data=memdata,aes(x=cellnum,y=memory,color=method))+
  scale_color_manual(values=colv2) +
  theme_minimal() + xlab(expression(paste('Number of cells (', "log"[10],'-scaled)')))+ ylab('Memory (GB)')
dev.off()
