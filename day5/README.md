## 1. 各种常用注释包简介
#### 1.1 TXDb
TxDb包，包含位置信息。包里的注释内容，可以通过包名来判断。
包名的格式是：TxDb.Species.Source.Build.Table
例如，TxDb.Hsapiens.UCSC.hg19.knownGene，代表：
+ 物种是 Homo sapiens
+ 来源于 UCSC genome browser
+ 基因组版本是 hg19 (their version of GRCh37)
+ 包含的注释内容是 knownGene table
此外还有，
TxDb.Dmelanogaster.UCSC.dm3.ensGene 
TxDb.Athaliana.BioMart.plantsmart22
等等

#### 1.2 EnsDb packages
EnsDb和TxDb类似，只不过它是基于Ensembl数据库的
例如：
+ EnsDb.Hsapiens.v79 
+ EnsDb.Mmusculus.v79 
+ EnsDb.Rnorvegicus.v79

#### 1.3 不常规使用示例
```r
select(TxDb.Hsapiens.UCSC.hg19.knownGene, c("1","10"),
       c("TXNAME","TXCHROM","TXSTART","TXEND"), "GENEID")
#> 'select()' returned 1:many mapping between keys and columns
#>   GENEID     TXNAME TXCHROM  TXSTART    TXEND
#> 1      1 uc002qsd.4   chr19 58858172 58864865
#> 2      1 uc002qsf.2   chr19 58859832 58874214
#> 3     10 uc003wyw.1    chr8 18248755 18258723
select(EnsDb.Hsapiens.v79, c("1", "10"),
       c("GENEID","GENENAME","SEQNAME","GENESEQSTART","GENESEQEND"), "ENTREZID")
#>   ENTREZID          GENEID GENENAME SEQNAME GENESEQSTART GENESEQEND
#> 1        1 ENSG00000121410     A1BG      19     58345178   58353499
#> 2       10 ENSG00000156006     NAT2       8     18391245   18401218
```

## 2. GRanges，Ranges
#### 2.1 GRanges
注释包比较常用的一个功能，是将注释包中的位置信息提取出来，添加到`GRanges`或者`GRangesList`对象中。 
以下示例查看注释包中基因的位置信息
```r
gns <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gns
#> GRanges object with 23056 ranges and 1 metadata column:
#>         seqnames              ranges strand |     gene_id
#>            <Rle>           <IRanges>  <Rle> | <character>
#>       1    chr19   58858172-58874214      - |           1
#>      10     chr8   18248755-18258723      + |          10
#>     100    chr20   43248163-43280376      - |         100
#>    1000    chr18   25530930-25757445      - |        1000
#>   10000     chr1 243651535-244006886      - |       10000
#>     ...      ...                 ...    ... .         ...
#>    9991     chr9 114979995-115095944      - |        9991
#>    9992    chr21   35736323-35743440      + |        9992
#>    9993    chr22   19023795-19109967      - |        9993
#>    9994     chr6   90539619-90584155      + |        9994
#>    9997    chr22   50961997-50964905      - |        9997
#>   -------
#>   seqinfo: 93 sequences (1 circular) from hg19 genome
```
#### 2.2 GRangesList
```r
## the genomic position of all transcripts by gene
txs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene)
txs
#> GRangesList object of length 23459:
#> $1 
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames            ranges strand |     tx_id     tx_name
#>          <Rle>         <IRanges>  <Rle> | <integer> <character>
#>   [1]    chr19 58858172-58864865      - |     70455  uc002qsd.4
#>   [2]    chr19 58859832-58874214      - |     70456  uc002qsf.2
#> 
#> $10 
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames            ranges strand | tx_id    tx_name
#>   [1]     chr8 18248755-18258723      + | 31944 uc003wyw.1
#> 
#> $100 
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames            ranges strand | tx_id    tx_name
#>   [1]    chr20 43248163-43280376      - | 72132 uc002xmj.3
#> 
#> ...
#> <23456 more elements>
#> -------
#> seqinfo: 93 sequences (1 circular) from hg19 genome
```
以List形式呈现，每个基因所对应的转录本，在基因组上的坐标

#### 2.3 其它用于注释包中选取信息的函数
+ 也可以使用，`transcripts()`, `genes()` 函数选取位置信息， 编码区可以使用`cds()`,`promoters()`,`exons()` 来选取
+ Positional information can be extracted for most of the above, grouped by a second element. For example, our `transcriptsBy ` call was all transcripts, grouped by gene

**`Ranges`** 对象功能十分强大，可以让我们基于基因组位置信息轻松选取出相关数据。
`GRanges`和`GRangesLists`对象，就类似于data.frame和List数据结构，我们可以使用`[`符号去选取子集。
示例
```r
txs[txs %over% gns[1:2,]]
#> GRangesList object of length 3:
#> $1 
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames            ranges strand |     tx_id     tx_name
#>          <Rle>         <IRanges>  <Rle> | <integer> <character>
#>   [1]    chr19 58858172-58864865      - |     70455  uc002qsd.4
#>   [2]    chr19 58859832-58874214      - |     70456  uc002qsf.2
#> 
#> $10 
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames            ranges strand | tx_id    tx_name
#>   [1]     chr8 18248755-18258723      + | 31944 uc003wyw.1
#> 
#> $162968 
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames            ranges strand | tx_id    tx_name
#>   [1]    chr19 58865723-58874214      - | 70457 uc002qsh.2
#>   [2]    chr19 58865723-58874214      - | 70458 uc002qsi.2
#> 
#> -------
#> seqinfo: 93 sequences (1 circular) from hg19 genome
```
 
## 3. TxDb包相关小练习
+ How many transcripts does PPARG have, according to UCSC?

+ Does Ensembl agree?

+ How many genes are between 2858473 and 3271812 on chr2 in the hg19 genome?