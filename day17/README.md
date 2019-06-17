[toc]

使用limma,edgeR,Glimma 进行完整的数据分析流程指南 可参见：
http://master.bioconductor.org/packages/release/workflows/html/RNAseq123.html

---
## 1. 学习目标
+ read in count data and format as a DGEList-object 读取count格式的数据，并格式化为DGE对象
+ annotate Entrez gene identifiers with gene information 注释基因
+ filter out lowly expressed genes 过滤低表达的基因
+ normalise gene expression values 标准化基因表达数据
+ unsupervised clustering of samples (standard and interactive plots) 对样本进行无监督聚类
+ linear modelling for comparisons of interest 对感兴趣的分组用线性模型进行比较
+ remove heteroscedascity
+ examine the number of differentially expressed genes 检查差异表达基因的数目
+ mean-difference plots (standard and interactive plots) 画图
+ heatmaps 绘制热图


#### 1. 安装RNAseq123包和所需的包，并下载样本数据
```r
suppressPackageStartupMessages({
  library(limma)
  library(Glimma)
  library(edgeR)
  library(Mus.musculus)
})
```

```r
BiocManager::install("RNAseq123")
## download sample data
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"  
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb")   
utils::untar("GSE63310_RAW.tar", exdir = ".")  
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))  
  R.utils::gunzip(i, overwrite=TRUE)  
```

#### 2. 读入文件
每份文件，都包含给定样本的count数；我们本次练习只包含 basal, LP 和 ML这三类样本
```r
## read data
read.delim(file.path(".", files[1]), nrow=5)

## Use readDEG read counts
x <- readDGE(file.path(".", files), columns=c(1,3))
class(x)
dim(x)
```
如果count的数据是存储在单个文件中，那么可以使用 `DEGList` 函数，将数据读入函数后转为 `DGEList` 对象

#### 3. 构建样本的分组信息







