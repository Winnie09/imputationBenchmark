### plot imputation
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
pd = readRDS('./plot/pd/11_imp_eval_4subfig.pd.rds')
pd4 = pd$pd4
pd4$platform = sapply(pd4[, 2], function(i) sub('_.*', '', i))

pval <- t(sapply(unique(pd4[,1]), function(method){
  tmp = pd4[pd4[,1]==method,]
  r = t.test(tmp[tmp[,4]=='fluidigm',3], tmp[tmp[,4]=='10x',3] )
  c(r$p.value, r$statistic, r$estimate)
}))
pval = cbind(pval, fdr = p.adjust(pval[,1], method = 'fdr'))
colnames(pval)[1] <- 'p.value'
rownames(pval) <- unique(pd4[,1])
pval = pval[order(pval[,'fdr']),]
write.csv(pval, './plot/pd/11_imp_eval_4subfig_compare_platform_ttest.csv')
#                 p.value          t  mean of x  mean of y          fdr
# scimpute   2.867148e-06 12.1942247 0.69171643 0.52795541 0.0000544758
# dca        1.212285e-05 10.4840166 0.74300573 0.57475754 0.0001151671
# scVI       1.951885e-05  9.5479742 0.74287931 0.57717273 0.0001236194
# pblr       7.563397e-05  8.7808483 0.64518994 0.45733837 0.0003592614
# mcimpute   1.813352e-04  7.2900321 0.65370246 0.50790390 0.0006890739
# deepimpute 4.442280e-04  8.8573545 0.58441574 0.44408284 0.0014067220
# autoimpute 5.890912e-04  9.2181666 0.31595212 0.07437382 0.0015989619
# magic      8.709780e-03  3.7608574 0.67615417 0.60076585 0.0206857285
# saucie     1.141392e-02  4.1882091 0.63863715 0.51687919 0.0240960623
# scscope    2.752200e-02 -3.3881401 0.08115045 0.33466438 0.0522918067
# viper      3.282514e-02  3.0409510 0.60189430 0.44936353 0.0566979617
# knnsmooth  4.605405e-02  2.7187143 0.67861236 0.55808838 0.0729189181
# drimpute   9.432221e-02  2.1009243 0.62516664 0.54005631 0.1378555321
# saverx     1.955402e-01 -1.4701257 0.54613587 0.58476234 0.2653759256
# saver      2.663916e-01 -1.2198874 0.56008749 0.58496730 0.3374294198
# alra       3.401775e-01  1.0640005 0.54237874 0.48021992 0.4039607223
# raw        3.620266e-01  1.0213457 0.48022447 0.41669407 0.4046179502
# baynorm    4.989340e-01 -0.7413000 0.48307031 0.52800903 0.5266525702
# screcover  7.967850e-01  0.2745928 0.54347520 0.52229593 0.7967849645
#     scimpute          dca         scVI         pblr     mcimpute   deepimpute 
# 0.0000544758 0.0001151671 0.0001236194 0.0003592614 0.0006890739 0.0014067220 
#   autoimpute        magic       saucie      scscope        viper    knnsmooth 
# 0.0015989619 0.0206857285 0.0240960623 0.0522918067 0.0566979617 0.0729189181 
#     drimpute       saverx        saver         alra          raw      baynorm 
# 0.1378555321 0.2653759256 0.3374294198 0.4039607223 0.4046179502 0.5266525702 
#    screcover 
# 0.7967849645 


v = sapply(1:nrow(pval), function(i) abs(pval[i,3] - pval[i,4])/max(pval[i,3], pval[i,4])) 
pval = cbind(pval, mean_improve = v)
rownames(pval)[pval[,'fdr'] < 0.05 & pval[,'mean_improve'] > 0.25]
pval1 = pval

###
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/rna_imputation/')
# setwd('/Users/wenpinhou/Dropbox/rna_imputation/')
source('./resource/function.R')
pd <- readRDS('./plot/pd/15_imp_diff_eval_4subfig.pd.rds')
pd4 = pd$pd4
pval <- t(sapply(unique(pd4[,1]), function(method){
  tmp = pd4[pd4[,1]==method,]
  r = t.test(tmp[tmp[,4]=='fluidigm',3], tmp[tmp[,4]=='10x',3] )
  c(r$p.value, r$statistic, r$estimate)
}))
pval <- cbind(pval, fdr = p.adjust(pval[,1], method = 'fdr'))
colnames(pval)[1] <- 'p.value'
rownames(pval) <- unique(pd4[,1])
pval <- pval[order(pval[,'fdr']),]
write.csv(pval, './plot/pd/15_imp_diff_eval_4subfig_compare_platform_ttest.csv')
#                 p.value           t    mean of x  mean of y          fdr
# saverx     1.076742e-11 -15.1995543  0.282923745 0.61476357 2.045809e-10
# saver      4.758672e-11 -13.8813952  0.270594982 0.60845484 4.051300e-10
# scscope    6.396789e-11 -14.0333788 -0.008467802 0.08546584 4.051300e-10
# autoimpute 9.542855e-06   7.7801623  0.058153012 0.01522261 4.532856e-05
# saucie     1.493266e-05   8.2923424  0.314259806 0.01208061 5.674411e-05
# magic      9.108459e-05  -6.0403747  0.451559391 0.70791327 2.884345e-04
# alra       1.720846e-04  -4.9003141  0.363517516 0.49925693 4.670869e-04
# drimpute   1.898145e-03  -3.6856349  0.338915290 0.42512074 4.508093e-03
# screcover  3.159237e-03  -3.4623356  0.272497723 0.37093746 6.669501e-03
# knnsmooth  5.939674e-03  -3.3388989  0.379860525 0.51014202 1.128538e-02
# scVI       6.577917e-03  -3.2019505  0.408740362 0.53782753 1.136186e-02
# dca        1.077789e-02  -2.8882978  0.447243152 0.53980219 1.575230e-02
# deepimpute 1.037501e-02  -2.8693781  0.191926514 0.24189973 1.575230e-02
# mcimpute   2.093427e-02  -2.5357856  0.340127590 0.39668255 2.841079e-02
# scimpute   2.256193e-02  -2.5397927  0.318774796 0.38988416 2.857844e-02
# baynorm    1.671119e-01  -1.4509156  0.234346533 0.26206243 1.984454e-01
# raw        3.393088e-01  -0.9858840  0.243886340 0.26225878 3.792275e-01
# pblr       7.240243e-01  -0.3586675  0.314051875 0.32379806 7.642479e-01
# viper      7.908992e-01   0.2694885  0.331459676 0.32528400 7.908992e-01

#       saverx        saver      scscope   autoimpute       saucie        magic 
# 2.045809e-10 4.051300e-10 4.051300e-10 4.532856e-05 5.674411e-05 2.884345e-04 
#         alra     drimpute    screcover    knnsmooth         scVI          dca 
# 4.670869e-04 4.508093e-03 6.669501e-03 1.128538e-02 1.136186e-02 1.575230e-02 
#   deepimpute     mcimpute     scimpute      baynorm          raw         pblr 
# 1.575230e-02 2.841079e-02 2.857844e-02 1.984454e-01 3.792275e-01 7.642479e-01 

v = sapply(1:nrow(pval), function(i) abs(pval[i,3] - pval[i,4])/max(pval[i,3], pval[i,4])) 
pval = cbind(pval, mean_improve = v)
rownames(pval)[pval[,'fdr'] < 0.05 & pval[,'mean_improve'] > 0.25]
pval2 = pval


###
rownames(pval1)[pval1[,'fdr']<0.05 & pval1[,'mean_improve'] > 0.25]
# [1] "pblr"       "autoimpute"
rownames(pval2)[pval2[,'fdr']<0.05 & pval2[,'mean_improve'] > 0.25]
# [1] "saverx"     "saver"      "scscope"    "autoimpute" "saucie"    
# [6] "magic"      "alra"       "screcover"  "knnsmooth"    
# autoimpute is better in 10x in both analysis


rownames(pval1)[pval1[,'fdr']<0.05 & pval1[,'mean_improve'] > 0.20]
# [1] "scimpute"   "dca"        "scVI"       "pblr"       "mcimpute"  
# [6] "deepimpute" "autoimpute"
rownames(pval2)[pval2[,'fdr']<0.05 & pval2[,'mean_improve'] > 0.20]
#  [1] "saverx"     "saver"      "scscope"    "autoimpute" "saucie"    
#  [6] "magic"      "alra"       "drimpute"   "screcover"  "knnsmooth" 
# [11] "scVI"       "deepimpute"

rownames(pval1)[pval1[,'fdr']<0.05 & pval1[,'mean_improve'] > 0.30]
# [1] "autoimpute"
rownames(pval2)[pval2[,'fdr']<0.05 & pval2[,'mean_improve'] > 0.30]
# [1] "saverx"     "saver"      "scscope"    "autoimpute" "saucie"    
# [6] "magic"     

rownames(pval1)[pval1[,'fdr']<0.05 & pval1[,'mean_improve'] > 0.3]
# [1] "autoimpute"
rownames(pval2)[pval2[,'fdr']<0.05 & pval2[,'mean_improve'] > 0.3]
# [1] "saverx"     "saver"      "scscope"    "autoimpute" "saucie"    
# [6] "magic"


##
rownames(pval1)[pval1[,'fdr']<0.05 & abs(pval1[,3]-pval1[,4]) > 0.15]
# [1] "scimpute"   "dca"        "scVI"       "pblr"       "autoimpute"
rownames(pval2)[pval2[,'fdr']<0.05 & abs(pval2[,3]-pval2[,4]) > 0.15]
# [1] "saverx" "saver"  "saucie" "magic" 
## no overlap


rownames(pval1)[pval1[,'fdr']<0.05 & abs(pval1[,3]-pval1[,4]) > 0.1]
# [1] "scimpute"   "dca"        "scVI"       "pblr"       "mcimpute"  
# [6] "deepimpute" "autoimpute" "saucie"
rownames(pval2)[pval2[,'fdr']<0.05 & abs(pval2[,3]-pval2[,4]) > 0.1]
# [1] "saverx"    "saver"     "saucie"    "magic"     "alra"      "knnsmooth"
# [7] "scVI" 
## scVI overlap and different direction

rownames(pval1)[pval1[,'fdr']<0.05 & abs(pval1[,3]-pval1[,4]) > 0.2]
# [1] "autoimpute"
rownames(pval2)[pval2[,'fdr']<0.05 & abs(pval2[,3]-pval2[,4]) > 0.2]
# [1] "saverx" "saver"  "saucie" "magic" 

