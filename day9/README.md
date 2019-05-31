[toc]

## 1. BSgenome packages
+ BSgenome包包含了给定物种/基因组版本的序列信息。我们可以通过使用 `available.genomes`函数，获取可用的BSgenome系列包
```r

library(BSgenome)
head(available.genomes())
#> [1] "BSgenome.Alyrata.JGI.v1"                 "BSgenome.Amellifera.BeeBase.assembly4"  
#> [3] "BSgenome.Amellifera.UCSC.apiMel2"        "BSgenome.Amellifera.UCSC.apiMel2.masked"
#> [5] "BSgenome.Athaliana.TAIR.04232008"        "BSgenome.Athaliana.TAIR.TAIR9"
```
+ 加载并检查BSgenome包
```r
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
```
+ 使用 **`getSeq()`** 函数获取序列数据，或者通过传入GRanges对象来获取某个区域的数据
```r
getSeq(Hsapiens, "chr1")
#>   249250621-letter "DNAString" instance
#> seq: NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
getSeq(Hsapiens, gns["5467",])
#>   A DNAStringSet instance of length 1
#>     width seq                                                                                       names               
#> [1] 85634 GCGGAGCGTGTGACGCTGCGGCCGCCGCGGACCTGGGGATTAA...ACTTTAAATAAATGGGAATTAAATATTTAAGAGCTGACTGGAA 5467
```
Biostrings包包含处理这些 `*StringSet` 对象的大部分代码—可以参阅Biostrings小片段和帮助页面以获得更多信息。
+ **BSgenome 练习** Get the sequences for all transcripts of the TP53 gene
````r

```

## 2. AnnotationHub
AnnotationHub是一个R包，也可以视作一个通道，我们可以在不安装相应注释包的情况下，通过它下载和查询多种不同的注释对象（annotation objects）
```r
library(AnnotationHub)
hub <- AnnotationHub()
#> snapshotDate(): 2018-06-27
hub
```

#### 2.1 查询AnnotationHub
想在AnnotationHub上找到“正确”的资源，就像使用谷歌一样，我们需要使用一个经过良好配置的查询检索方法，才能获得自己想要的信息。
在`AnnotationHub`中进行高效的查询，需要注意根据以下规则分类检索：
+ Data provider
+ Data class
+ Species
+ Data source
```r
names(mcols(hub))
#>  [1] "title"              "dataprovider"       "species"            "taxonomyid"         "genome"            
#>  [6] "description"        "coordinate_1_based" "maintainer"         "rdatadateadded"     "preparerclass"     
#> [11] "tags"               "rdataclass"         "rdatapath"          "sourceurl"          "sourcetype"
```

#### 2.2 AnnotationHub Data providers, Data classes, Species， Data sources
```r
## AnnotationHub Data providers
unique(hub$dataprovider)

## AnnotationHub Data classes
unique(hub$rdataclass)

## AnnotationHub Species
head(unique(hub$species))
length(unique(hub$species))

## Data sources
unique(hub$sourcetype)
```
#### 2.3 AnnotationHub query
```r
## AnnotationHub query
qry <- query(hub, c("granges","homo sapiens","ensembl"))
qry
qry$sourceurl

## Selecting AnnotationHub resource
whatIwant <- qry[["AH50377"]]
GRCh38TxDb <- makeTxDbFromGRanges(whatIwant)
GRCh38TxDb
```

#### 2.4 AnnotationHub exercises
+ How many resources are on AnnotationHub for Atlantic salmon (Salmo salar)?
+ Get the most recent Ensembl build for domesticated dog (Canis familiaris) and make a TxDb


#### 2.5 AnnotationHub使用的报错解决
+  重要参考网址：http://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/TroubleshootingTheCache.html