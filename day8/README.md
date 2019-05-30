[toc]
## 学习章节
https://bioconductor.github.io/BiocWorkshops/introduction-to-bioconductor-annotation-resources.html#organism.dplyr-package

## 1. OrganismDb.dplyr 包
> + 简介：The Organism.dplyr creates an on disk sqlite database to hold data of an organism combined from an ‘org’ package (e.g., org.Hs.eg.db) and a genome coordinate functionality of the ‘TxDb’ package
> + It aims to provide an integrated presentation of identifiers and genomic coordinates. And a src_organism object is created to point to the database.
> + The *src_organism* object is created as an extension of *src_sql* and *src_sqlite* from `dplyr`, which inherited all `dplyr` methods. It also implements the `select() ` interface from `AnnotationDbi` and genomic coordinates extractors from `GenomicFeatures`


+ OrganismDb.dplyr包将TxDb和Org.Db包中的相关数据整合到了本地
+ 可以使用来自`Org.*` 和 `TxDb.*`的函数提取注释包中的数据
  + `keytypes()`, `select()`, ...
  + `exons()`, `promoters()`, ...
+ 可以使用`dplyr`的函数来显示过滤 **`TxDb`**和**`Org.Db`**的整合注释
```r
library(Organism.dplyr)
src <- src_organism(dbpath = hg38light())
src
## src:  sqlite 3.22.0 [D:\software\R-3.5.3\library\Organism.dplyr\extdata\light.hg38.knownGene.sqlite]
## tbls: id, id_accession, id_go, id_go_all, id_omim_pm,
##  id_protein, id_transcript, ranges_cds, ranges_exon,
##  ranges_gene, ranges_tx
```
+ **`解释-hg38light()`**
> These functions are primarily for illustrating functionality.hg38light() and mm10light() provide access to
> trimmed-down versions of Organism.dplyr data based derived from the TxDb.Hsapiens.UCSC.hg38.knownGene and 
> TxDb.Mmusculus.UCSC.mm10.ensGene data bases

+ **`src_organism()`**
> Create a sqlite database from TxDb and corresponding Org packages


#### 1.0 基础用法
+ 建议使用 **`browseVignettes("Organism.dplyr")`** 来查看这个包的具体使用方法
```r
## Look at all available tables
src_tbls(src)
##  [1] "id_accession"  "id_transcript" "id"            "id_omim_pm"   
##  [5] "id_protein"    "id_go"         "id_go_all"     "ranges_gene"  
##  [9] "ranges_tx"     "ranges_exon"   "ranges_cds"

## Look at data from one specific table
tbl(src, "id")
## # Source:   table<id> [?? x 6]
## # Database: sqlite 3.22.0 []
##    entrez map      ensembl         symbol genename               alias   
##    <chr>  <chr>    <chr>           <chr>  <chr>                  <chr>   
##  1 1      19q13.4  ENSG00000121410 A1BG   alpha-1-B glycoprotein A1B     
##  2 1      19q13.4  ENSG00000121410 A1BG   alpha-1-B glycoprotein ABG     

## Look at fields of one table.
colnames(tbl(src, "id"))
## [1] "entrez"   "map"      "ensembl"  "symbol"   "genename" "alias"
```

#### 1.1 获取TxDb对象中的启动子
```r
options(width = 120)
promoters(src)
```

#### 1.2 Extract a table from the underlying database
```r
tbl(src, "id")
#> # Source:   table<id> [?? x 6]
#> # Database: sqlite 3.22.0 []
#>    entrez map      ensembl         symbol genename               alias   
#>    <chr>  <chr>    <chr>           <chr>  <chr>                  <chr>   
#>  1 1      19q13.4  ENSG00000121410 A1BG   alpha-1-B glycoprotein A1B     
#>  2 1      19q13.4  ENSG00000121410 A1BG   alpha-1-B glycoprotein ABG     
#>  3 1      19q13.4  ENSG00000121410 A1BG   alpha-1-B glycoprotein GAB     
#>  4 1      19q13.4  ENSG00000121410 A1BG   alpha-1-B glycoprotein HYST2477
```
#### 1.3 执行一个复杂查询 Make a complex query between tables in the underlying database
```r
inner_join(tbl(src, "id"), tbl(src, "ranges_gene")) %>%
            filter(symbol %in% c("ADA", "NAT2")) %>%
            dplyr::select(gene_chrom, gene_start, gene_end,
            gene_strand, symbol, alias, map)
```

#### 1.4 Organism.dplyr 练习
+ How many supported organisms are implemented in Organism.dplyr?
```r
browseVignettes("Organism.dplyr")
supportedOrganisms()
## 21
```
+ Display the ensembl Id and genename description for symbol “NAT2”
```r
# look at all avaliable tables
src_tbls(src)
# get ensembl Id and genename description
 tbl(src, "id") %>%
  filter(symbol == "NAT2") %>%
  dplyr::select( ensembl,genename)
```
+ Show all the alias for “NAT2” in the database
```r
tbl(src, "id") %>%
  filter(symbol == "NAT2") %>%
  dplyr::select(alias)
```

+ Display Gene ontology (GO) information for gene symbol “NAT2”
```r
inner_join(tbl(src, "id"), tbl(src, "id_go")) %>%
    filter(symbol == "NAT2") %>%
    dplyr::select(entrez, ensembl, symbol, go, evidence, ontology)
```

+ Gene transcript counts per gene symbol
```r
txcount <- inner_join(tbl(src, "id"), tbl(src, "ranges_tx")) %>%
    dplyr::select(symbol, tx_id) %>%
    group_by(symbol) %>%
    summarize(count = n()) %>%
    dplyr::select(symbol, count) %>%
    arrange(desc(count)) %>%
    collect(n=Inf)
txcount
## plot
library(ggplot2)
ggplot(txcount, aes(x = symbol, y = count)) + 
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Transcript count") +
    labs(x = "Symbol") +
    labs(y = "Count")

```