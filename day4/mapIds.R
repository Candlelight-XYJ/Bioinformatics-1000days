library(hugene20sttranscriptcluster.db)

ids <- c('16737401','16657436' ,'16678303')
## List
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID", multiVals = "list")
## CharacterList
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID", multiVals = "CharacterList")
## filter
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID", multiVals = "filter")
## asNA
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID", multiVals = "asNA")

##################
#### practice ####
##################

## What gene symbol corresponds to Entrez Gene ID 1000?
entrezid="1000"
columns(hugene20sttranscriptcluster.db)
select(hugene20sttranscriptcluster.db,entrezid,"SYMBOL",keytype = "ENTREZID")
# 'select()' returned 1:1 mapping between keys and columns
# ENTREZID SYMBOL
# 1     1000   CDH2

## What is the Ensembl Gene ID for PPARG?
select(hugene20sttranscriptcluster.db,"PPARG","ENSEMBL",keytype = "SYMBOL")
# 'select()' returned 1:1 mapping between keys and columns
# SYMBOL         ENSEMBL
# 1  PPARG ENSG00000132170  

## What is the UniProt ID for GAPDH?
select(hugene20sttranscriptcluster.db,"GAPDH","UNIPROT",keytype = "SYMBOL")
# 'select()' returned 1:many mapping between keys and columns
# SYMBOL UNIPROT
# 1  GAPDH  P04406
# 2  GAPDH  V9HVZ4

## How many of the probesets from the ExpressionSet (eset) we loaded map to a single gene? 

load("E:\\GitHub\\Bioinformatics-1000days\\day3\\eset.Rdata")  
#genes <- as.character(pData(featureData(eset))$SYMBOL)
#res <- mapIds(hugene20sttranscriptcluster.db,genes,"PROBEID",keytype = "SYMBOL",multiVals = "CharacterList")
ids <- featureNames(eset)
res <- mapIds(hugene20sttranscriptcluster.db,ids,"SYMBOL","PROBEID",multiVals = "filter")
length(res)
#[1] 31792

## How many don’t map to a gene at all?
res <- mapIds(hugene20sttranscriptcluster.db,ids,"SYMBOL","PROBEID",multiVals = "CharacterList")
df <- data.frame(id=names(res),length=lengths(res))
num <- nrow(df[which(df$length>1),])
num
#[1] 1760




