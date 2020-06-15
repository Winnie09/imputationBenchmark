### plot imputation
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
pd = readRDS('./plot/pd/11_imp_eval_4subfig.pd.rds')
pd4 = pd$pd4
pd4$platform = sapply(pd4[, 2], function(i) sub('_.*', '', i))

pval <- sapply(unique(pd4[,1]), function(method){
  tmp = pd4[pd4[,1]==method,]
  t.test(tmp[tmp[,4]=='fluidigm',3], tmp[tmp[,4]=='10x',3] )$p.value
})
fdr = p.adjust(pval, method = 'fdr')
names(fdr) <- unique(pd4[,1])
fdr1 = fdr
# > sort(fdr1)
#     scimpute          dca         scVI         pblr     mcimpute   deepimpute 
# 0.0000544758 0.0001151671 0.0001236194 0.0003592614 0.0006890739 0.0014067220 
#   autoimpute        magic       saucie      scscope        viper    knnsmooth 
# 0.0015989619 0.0206857285 0.0240960623 0.0522918067 0.0566979617 0.0729189181 
#     drimpute       saverx        saver         alra          raw      baynorm 
# 0.1378555321 0.2653759256 0.3374294198 0.4039607223 0.4046179502 0.5266525702 
#    screcover 
# 0.7967849645 


setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
pd = readRDS('./plot/pd/15_imp_diff_eval_4subfig.pd.rds')
pd4 = pd$pd4
pval <- sapply(unique(pd4[,1]), function(method){
  tmp = pd4[pd4[,1]==method,]
  t.test(tmp[tmp[,4]=='fluidigm',3], tmp[tmp[,4]=='10x',3] )$p.value
})
fdr = p.adjust(pval, method = 'fdr')
names(fdr) <- unique(pd4[,1])
fdr2 = fdr
# > sort(fdr2)
#       saverx        saver      scscope   autoimpute       saucie        magic 
# 2.045809e-10 4.051300e-10 4.051300e-10 4.532856e-05 5.674411e-05 2.884345e-04 
#         alra     drimpute    screcover    knnsmooth         scVI          dca 
# 4.670869e-04 4.508093e-03 6.669501e-03 1.128538e-02 1.136186e-02 1.575230e-02 
#   deepimpute     mcimpute     scimpute      baynorm          raw         pblr 
# 1.575230e-02 2.841079e-02 2.857844e-02 1.984454e-01 3.792275e-01 7.642479e-01 
