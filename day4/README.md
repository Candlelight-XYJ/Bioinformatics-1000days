## 学习章节

## 1. The mapIdsfunction
**`mapIds()`** 函数与 **`select()`** 函数功能相似，但它可以对数据中出现的重复元素进行清理
`mapIds`与`select`函数在参数上的不同之处在于：
+ **`columns`** 参数，只能选择一个column
+ **`keytype`** 参数必须指定
+ **`mapIds`** 包含`multiVals`函数用于清理重复元素
```r
mapIds(hugene20sttranscriptcluster.db,ids,"SYMBOL","PROBEID")
```
#### 1.1 multiVals的使用
multiVals的默认参数是`first`（只选取重复元素中的第一个）。还有其他可选的选项：`list`, `CharacterList`,`filter`,`asNA`,也可以选择用户自定义的选项
+ **first:** This value means that when there are multiple matches only the 1st thing that comes back will be returned. This is the default behavior
+ **list:** This will just returns a list object to the end user
+ **filter:** This will remove all elements that contain multiple matches and will therefore return a shorter vector than what came in whenever some of the keys match more than one value
+ **asNA:** This will return an NA value whenever there are multiple matches
+ **CharacterList:** This just returns a SimpleCharacterList object
```r
## 选择list,就是以list的方式显示map结果
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID", multiVals = "list")
#> 'select()' returned 1:many mapping between keys and columns
#> $`16737401`
#> [1] "TRAF6"
#> 
#> $`16657436`
#> [1] "DDX11L1"      "LOC102725121" "DDX11L2"      "DDX11L9"     
#> [5] "DDX11L10"     "DDX11L5"      "DDX11L16"    
#> 
#> $`16678303`
#> [1] "ARF1"

## 选择CharacterList
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID", multiVals = "CharacterList")
'select()' returned 1:many mapping between keys and columns
CharacterList of length 3
[["16737401"]] TRAF6
[["16657436"]] DDX11L1 LOC102725121 DDX11L2 DDX11L9 DDX11L10 DDX11L5 DDX11L16
[["16678303"]] ARF1

## filter
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID", multiVals = "filter")
'select()' returned 1:many mapping between keys and columns
16737401 16678303 
 "TRAF6"   "ARF1"
## asNA
mapIds(hugene20sttranscriptcluster.db, ids, "SYMBOL", "PROBEID", multiVals = "asNA")
#> 'select()' returned 1:many mapping between keys and columns
#> 16737401 16657436 16678303 
#>  "TRAF6"       NA   "ARF1"
```
#### 1.2 小练习
+ What gene symbol corresponds to Entrez Gene ID 1000?
```r
entrezid="1000"
columns(hugene20sttranscriptcluster.db)
select(hugene20sttranscriptcluster.db,entrezid,"SYMBOL",keytype = "ENTREZID")
# 'select()' returned 1:1 mapping between keys and columns
# ENTREZID SYMBOL
# 1 1000 CDH2
```
+ What is the Ensembl Gene ID for PPARG?
```r
select(hugene20sttranscriptcluster.db,"PPARG","ENSEMBL",keytype = "SYMBOL")
# 'select()' returned 1:1 mapping between keys and columns
# SYMBOL ENSEMBL
# 1 PPARG ENSG00000132170  
```
+ What is the UniProt ID for GAPDH?
```r
select(hugene20sttranscriptcluster.db,"GAPDH","UNIPROT",keytype = "SYMBOL")
# 'select()' returned 1:many mapping between keys and columns
# SYMBOL UNIPROT
# 1 GAPDH P04406
# 2 GAPDH V9HVZ4
```
+ How many of the probesets from the ExpressionSet (eset) we loaded map to a single gene? 
```r
load("E:\\GitHub\\Bioinformatics-1000days\\day3\\eset.Rdata")  
#genes <- as.character(pData(featureData(eset))$SYMBOL)
#res <- mapIds(hugene20sttranscriptcluster.db,genes,"PROBEID",keytype = "SYMBOL",multiVals = "CharacterList")
ids <- featureNames(eset)
res <- mapIds(hugene20sttranscriptcluster.db,ids,"SYMBOL","PROBEID",multiVals = "filter")
length(res)
#[1] 31792
```
+ How many don’t map to a gene at all?
```r
## How many don’t map to a gene at all?
res <- mapIds(hugene20sttranscriptcluster.db,ids,"SYMBOL","PROBEID",multiVals = "CharacterList")
df <- data.frame(id=names(res),length=lengths(res))
num <- nrow(df[which(df$length>1),])
num
#[1] 1760
```