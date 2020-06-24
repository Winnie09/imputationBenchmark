setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
pd = readRDS('./plot/pd/11_imp_eval_4subfig.pd.rds')
coldf = readRDS('./resource/method_color.rds')
pd4 = pd[['pd4']]
cellnum = as.numeric(gsub('\\)','',sapply(pd4[,2],function(i) sub('.*\\(','',i))))
m = data.frame(pd4, cellnum = cellnum, stringsAsFactors = FALSE)

m = sapply(unique(pd4[,2]), function(i) pd4[which(pd4[,2]==i), 'cor'])
rownames(m) = pd4[which(pd4[,2]==unique(pd4[,2])[1]),1]
m1 = m[, 1:5]
m2 = m[, 6:10]

v1 = as.numeric(gsub('\\)', '', sapply(colnames(m1), function(i) sub('.*\\(', '', i))))
v2 = as.numeric(gsub('\\)', '', sapply(colnames(m2), function(i) sub('.*\\(', '', i))))

cor1 = apply(m1, 1, cor, v1, method = 'spearman')
cor2 = apply(m2, 1, cor, v2, method = 'spearman')

pd = readRDS('./plot/pd/15_imp_diff_eval_4subfig.pd.rds')
pd4 = pd[['pd4']]
# pd4[,'Method'] = factor(coldf[match(pd4[,'Method'],coldf[,'shortName']),'fullName'],levels=coldf[match(levels(pd4$Method),coldf[,'shortName']),'fullName'])

tmp1 = pd4[pd4[,4] == 'fluidigm', ]
tmp2 = pd4[pd4[,4] == '10x', ]
m1 = sapply(unique(tmp1[,2]), function(i) tmp1[which(tmp1[,2]==i), 3])
rownames(m1) = tmp1[which(tmp1[,2]==unique(tmp1[,2])[1]),1]
colnames(m1) = unique(tmp1[,2])
m2 = sapply(unique(tmp2[,2]), function(i) tmp2[which(tmp2[,2]==i), 3])
rownames(m2) = tmp2[which(tmp2[,2]==unique(tmp2[,2])[1]),1]
colnames(m2) = unique(tmp2[,2])

num1_1 = as.numeric(sapply(colnames(m1), function(i){
  sub('.*\\(', '', sub("\\)_.*", '', i))
}))
num1_2 = as.numeric(sapply(colnames(m1), function(i){
  sub('\\)', '', sub(".*\\(", '', i))
}))
num2_1 = as.numeric(sapply(colnames(m2), function(i){
  sub('.*\\(', '', sub("\\)_.*", '', i))
}))
num2_2 = as.numeric(sapply(colnames(m2), function(i){
  sub('\\)', '', sub(".*\\(", '', i))
}))

for (i in 1:length(num1_1)){
  if (num1_1[i] > num1_2[i]){
    tmp = num1_1[i]
    num1_1[i]  = num1_2[i]
    num1_2[i] = tmp
  }
}
for (i in 1:length(num1_1)){
  if (num2_1[i] > num2_2[i]){
    tmp = num2_1[i]
    num2_1[i] = num2_2[i]
    num2_2[i]  = tmp
  }
}
  
v11 = apply(m1, 1, cor, num1_1, method = 'spearman')
v12 = apply(m1, 1, cor, num1_2, method = 'spearman')
v13 = apply(m1, 1, cor, abs(num1_1 - num1_2), method = 'spearman')
v21 = apply(m2, 1, cor, num2_1, method = 'spearman')
v22 = apply(m2, 1, cor, num2_2, method = 'spearman')
v23 = apply(m2, 1, cor, abs(num2_1 - num2_2), method = 'spearman')

par(mfrow=c(2,3))
hist(v11, col = 'grey', main = 'fluidigm LFC & min.cellnum')
hist(v12, col = 'grey', main = 'fluidigm LFC & max.cellnum')
hist(v13, col = 'grey', main = 'fluidigm LFC & cellnum.diff')
hist(v21, col = 'grey', main = '10x LFC & min.cellnum')
hist(v22, col = 'grey', main = '10x LFC & max.cellnum')
hist(v23, col = 'grey', main = '10x LFC & cellnum.diff')

mtd = names(cor1)
df = data.frame(fluidigm_sc_score_and_cellnum_SCC = cor1,
                 umi_sc_score_and_cellnum_SCC = cor2[mtd],
                 fluidigm.LFC_score_and_min.cellnum_SCC = v11[mtd], 
                 fludigm.LFC_score_and_max.cellnum_SCC = v12[mtd], 
                 fluidigm.LFC_score_and_cellnum.diff_SCC = v13[mtd],
                 umi.LFC_score_and_min.cellnum_SCC = v21[mtd],
                 umi.LFC_score_and_max.cellnum_SCC = v22[mtd],
                 umi.LFC_score_and_cellnum.diff_SCC = v23[mtd], stringsAsFactors = FALSE)
write.csv(df, './plot/plot/SCC_between_cellnumber_bulkCorrelation.csv')

write.csv(df[, c(1,2,3,6)], './plot/plot/Fig2_SCC_between_cellnumber_bulkCorrelation.csv')
