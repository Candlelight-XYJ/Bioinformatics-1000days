## 学习章节
https://bioconductor.github.io/BiocWorkshops/r-and-bioconductor-for-everyone-an-introduction.html#introduction-to-bioconductor

[toc]
## 1. Bioconductor的一些补充小用法
使用valid()查看包安装的版本情况
```r
> BiocManager::valid()

* sessionInfo()

R version 3.5.3 (2019-03-11)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.936  LC_CTYPE=Chinese (Simplified)_China.936   
[3] LC_MONETARY=Chinese (Simplified)_China.936 LC_NUMERIC=C                              
[5] LC_TIME=Chinese (Simplified)_China.936    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] BiocManager_1.30.4 compiler_3.5.3     tools_3.5.3        yaml_2.2.0        

Bioconductor version '3.8'

  * 45 packages out-of-date
  * 0 packages too new

create a valid installation with

  BiocManager::install(c(
    "AnnotationFilter", "backports", "biomaRt", "Biostrings", "broom", "checkmate", "devtools",
    "digest", "dplyr", "DT", "ensembldb", "fs", "GenomicAlignments", "GenomicFeatures",
    "ggplot2", "ggrepel", "gridGraphics", "httpuv", "igraph", "knitr", "pillar", "processx",
    "progress", "ProtGenerics", "rcmdcheck", "RcppArmadillo", "remotes", "reprex", "rGREAT",
    "rlang", "robustbase", "rtracklayer", "rvest", "segmented", "shiny", "shinyFiles", "sys",
    "testthat", "tinytex", "usethis", "vegan", "WGCNA", "xfun", "xtable", "zip"
  ), update = TRUE, ask = FALSE)

more details: BiocManager::valid()$too_new, BiocManager::valid()$out_of_date

Warning message:
45 packages out-of-date; 0 packages too new 

```
avaliable() 用于搜索相关的包
```r
## 例如这里输入TxDb.Hsapiens，它会自动匹配相关的包
BiocManager::available("TxDb.Hsapiens")
#> [1] "TxDb.Hsapiens.BioMart.igis"                 
#> [2] "TxDb.Hsapiens.UCSC.hg18.knownGene"          
#> [3] "TxDb.Hsapiens.UCSC.hg19.knownGene"          
#> [4] "TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts"
#> [5] "TxDb.Hsapiens.UCSC.hg38.knownGene"
```
通过网站https://support.bioconductor.org/ 来查询解决和包相关的问题

## 2. Working with Genomic Ranges
+ Importing data

`rtracklayer`包提供函数`import()` 帮助用户读取基因组格式的文件(例如：BED,GTF,VCF,FASTA) 并封装成Bioconductor下的对象。`GenomicRanges`包提供了多种函数，来在基因组坐标系下操纵各种数据。
#### 2.1 importing data 导入bed文件
```r
library(rtracklayer)
## 使用file.choose() 来选择文件
fname <- file.choose()   # CpGislands.Hsapiens.hg38.UCSC.bed
fname
file.exists(fname)
## 使用import()函数导入bed文件。导入后，将以GenomicRanges的对象描述这份CpG岛数据。
cpg <- import(fname)
cpg
```
**`注意`**  BED格式的文件，它们的坐标体系是0-based的，并且intervals是半开区间(start在范围内，end在范围后一个坐标里)
而Bioconductor得坐标是1-based，并且是闭区间(start和end坐标都包含在范围内)，因此使用import()函数导入数据得时候，它会自动将bed文件坐标转换为Bioconductor文件对象的坐标。

#### 2.2 Working with genomic ranges

```r
## 使用函数 keepStandardChromosomes 保留标准染色体，标准染色体指的是chr1-22和x,y染色体
## 很多时候我们获取的数据里 可能染色体不止chr1-22与x,y ； 
## 可能还包含其它染色体，例如：chr22_KI270738v1_random这样的染色体。
## 通过keepStandardChromosomes 就可以去除这些其它的染色体，只保留标准染色体
cpg <- keepStandardChromosomes(cpg, pruning.mode = "coarse")
cpg
```
GenomicRanges对象包含两个部分
+ 必需的：seqnames（染色体号），start（起始位点）,end(终止位点)，strand (正链或负链)
+ 非必需的：另外的元素，例如本例中的name。

必需的元素内容，可以使用函数`start()`, `end()`, `width()`获取。非必需的元素，使用`$`符号获取内容。

```r
head( start(cpg) )
#> [1] 155188537   2226774  36306230  47708823  53737730 144179072
head( cpg$name )

## 使用subset()函数获取染色体1号和2号上的cpg岛
subset(cpg, seqnames %in% c("chr1", "chr2"))
```
#### 2.3 Genomic annotations
```r
## 导入注释数据
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
## 查看注释数据中的所有转录本
tx <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
tx
## 为了演示方便，目前只保留标准染色体
tx <- keepStandardChromosomes(tx, pruning.mode="coarse")
tx
```

#### 2.4 Overlaps
可以使用`findOverlaps()`,`nearest()`,`precede()`和`follow()`函数，来完成overlaps
```r
## count the number of CpG islands that overlap each transcript
olaps <- countOverlaps(tx, cpg)
length(olaps)     # 1 count per transcript
#> [1] 182435
table(olaps)
## 将计算的overlaps加入GR对象 
tx$cpgOverlaps <- countOverlaps(tx, cpg)
tx
## subsetting the GRanges objects for transcripts satisfying particular conditions, in a coordinated fashion ## where the software does all the book-keeping to makes sure the correct ranges are selected.
subset(tx, cpgOverlaps > 10)

```