setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation')
library(RColorBrewer)
library(ggplot2)
mtd = read.table('./resource/impute_method.txt')
fn = c('no_imp','ALRA','AutoImpute','bayNorm','DCA','DeepImpute','DrImpute','kNN-smoothing','MAGIC','mcImpute','PBLR','SAUCIE','SAVER','SAVERX','scImpute','scRecover','scScope','scVI','VIPER')
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(fn))
v = getPalette(colourCount)
v =c(setdiff(v,"#999999"),"#67001F")
df = data.frame(shortName=mtd,fullName=fn,col=v, stringsAsFactors = F)
colnames(df) = c('shortName','fullName','color')
df[df$shortName=='saucie','color'] = 'darkolivegreen'
df[df$shortName=='saverx','color'] = 'deeppink3'
# pd = cbind(df,den=1:length(v))
# ggplot(pd) + geom_histogram(aes(x=shortName), fill=pd$color,stat='count') +
# theme(legend.position="right")+ theme(axis.text.x = element_text(angle=45,hjust = 1))
saveRDS(df,'./resource/method_color.rds')


mtd = read.table('./resource/impute_method_latent.txt',stringsAsFactors = F)
fn = c('no_imp','ALRA','AutoImpute','bayNorm','DCA','DeepImpute','DrImpute','kNN-smoothing','MAGIC','mcImpute','PBLR','SAUCIE','SAVER','SAVERX','scImpute','scRecover','scScope','scVI','VIPER','SAUCIE_latent','scScope_latent','scVI_latent')
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(fn))-3
v = getPalette(colourCount)
v =c(setdiff(v,"#999999"),"#67001F")
v = c(v,'yellow4','navy','purple')
df = data.frame(shortName=mtd,fullName=fn,col=v, stringsAsFactors = F)
colnames(df) = c('shortName','fullName','color')
df[df$shortName=='saucie','color'] = 'darkolivegreen'
df[df$shortName=='saverx','color'] = 'deeppink3'
# pd = cbind(df,den=1:length(v))
# pd$shortName = factor(pd$shortName,levels=pd$shortName)
# ggplot(pd) + geom_histogram(aes(x=shortName), fill=pd$color,stat='count') +
#   theme(legend.position="right")+ theme(axis.text.x = element_text(angle=45,hjust = 1))
saveRDS(df,'./resource/method_latent_color.rds')




