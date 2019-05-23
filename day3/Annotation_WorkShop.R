setwd("E:/GitHub/Bioinformatics-1000days/day3/")
load("eset.Rdata")
eset

## use exprs() to check the expression data
head(exprs(eset))
## use pData() and phenoData() to check the pheno data
head(pData(phenoData(eset)))
## use pData() and featureData to check the featureData
## 例如基因的id,symbol之类的信息
head(pData(featureData(eset)))

########################################
##  Interacting with AnnoDb packages  ##
########################################

library(hugene20sttranscriptcluster.db)
set.seed(12345)
ids <- featureNames(eset)[sample(1:25000, 5)]
ids

select(hugene20sttranscriptcluster.db, ids, "SYMBOL")

## get keytypes or columns which are available for a given annotation package

keytypes(hugene20sttranscriptcluster.db)
columns(hugene20sttranscriptcluster.db)

## There is one issue with select however.
ids <- c('16737401','16657436' ,'16678303')
select(hugene20sttranscriptcluster.db, ids, c("SYMBOL","MAP"))


