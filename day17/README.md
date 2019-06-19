[toc]

使用limma,edgeR,Glimma 进行完整的数据分析流程指南 可参见：
http://master.bioconductor.org/packages/release/workflows/html/RNAseq123.html

---

## 1. 数据准备

####  学习目标
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


#### 1.1 安装RNAseq123包和所需的包，并下载样本数据
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

#### 1.2 读入文件
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

#### 1.3 构建样本的分组信息
对于下游分析而言，每个样本的分组实验信息需要添加上去。
包括： 细胞类型(本次练习中是 basal,LP,ML)，基因型(wild-type野生型，knock-out 敲除型)，表型(disease,status,sex,age)，样本处理情况(drug,control)，实验的批次信息等等
我们的DGEList对象包括：一个存储了所有细胞类型(`group`)和批次(`lane`)的 `samples` 的数据框。 注意，通过使用 `x$sanples` 提取数据时，每个样本的文库的大小会自动计算，并且 标准化因子会设置为1 。

+ 为了简化后续的运算，我们首先将GEO 样本ID中的 `GSM前缀` 去除
```r
## Organising sample information 设定sample的名字
# 提取sample的名字
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames

## 设置样本信息
colnames(x) <- samplenames
# 分组
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group ＃ 为每个样本添加细胞类型信息
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2))) 
x$samples$lane <- lane # 添加lane的信息
x$samples  # 此时sample信息中就会多出分组和细胞类型两列
```

#### 1.4 基因注释
**`DGEList`** 对象中包含名为 `genes` 的二级数据框，它主要用来存储和count矩阵相对应的行的基因相关信息。基因的信息可以使用物种包(例如小鼠的 `Mus.musculus` 包，人类的 `Homo.sapiens` 包) ，或者 `biomaRt` 包 来完成填充注释。


+ 本例中，我们使用提取物种包 `Mus.musculus` 中的注释基因 **gene symbols** 和染色体信息 来为count矩阵做注释（我们count矩阵原来是按照entrez id 来对基因进行表示的，这次我们通过与entrez id 对应的 gene symbol和染色体信息，来完善矩阵的基因注释）

```r
## using the Mus.musculus package 
## to retrieve associated gene symbols and chromosome information
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")   # 提取symbol和chr, 以entrezid作为map id的源头
head(genes) 
## !duplicated replicate genes
genes <- genes[!duplicated(genes$ENTREZID),]
## add genes 将提取好的基因信息添加入我们的DGEList中genes这个二级数据框
x$genes <- genes
x
```
**`需要注意的是`**
在本例中，注释和数据对象中的基因顺序是相同的。如果由于缺失和/或重新排列的基因id而不是这种情况，则可以使用 **`match`** 函数对基因进行正确排序。然后，将基因注释的数据框架添加到数据对象中，并整齐地打包在DGEList对象中，DGEList对象包含`原始计数数据`以及相关的`样本信息`和`基因注释`。

+ 这个时候就可以看到我们的DGEList对象x中有 `三` 个二级数据框了：`samples`， `counts`， `genes`

```r
> x
An object of class "DGEList"
$samples
                                 files group lib.size norm.factors lane
5_10_6_5_11 ./GSM1545535_10_6_5_11.txt    LP 32863052            1 L004
6_9_6_5_11   ./GSM1545536_9_6_5_11.txt    ML 35335491            1 L004
8_purep53     ./GSM1545538_purep53.txt Basal 57160817            1 L004
9_JMS8-2       ./GSM1545539_JMS8-2.txt Basal 51368625            1 L006
0_JMS8-3       ./GSM1545540_JMS8-3.txt    ML 75795034            1 L006
1_JMS8-4       ./GSM1545541_JMS8-4.txt    LP 60517657            1 L006
2_JMS8-5       ./GSM1545542_JMS8-5.txt Basal 55086324            1 L006
4_JMS9-P7c   ./GSM1545544_JMS9-P7c.txt    ML 21311068            1 L008
5_JMS9-P8c   ./GSM1545545_JMS9-P8c.txt    LP 19958838            1 L008

$counts
           Samples
Tags        5_10_6_5_11 6_9_6_5_11 8_purep53 9_JMS8-2 0_JMS8-3 1_JMS8-4 2_JMS8-5 4_JMS9-P7c
  497097              1          2       342      526        3        3      535          2
  100503874           0          0         5        6        0        0        5          0
  100038431           0          0         0        0        0        0        1          0
  19888               0          1         0        0       17        2        0          1
  20671               1          1        76       40       33       14       98         18
           Samples
Tags        5_JMS9-P8c
  497097             0
  100503874          0
  100038431          0
  19888              0
  20671              8
27174 more rows ...

$genes
   ENTREZID  SYMBOL TXCHROM
1    497097    Xkr4    chr1
2 100503874 Gm19938    <NA>
3 100038431 Gm10568    <NA>
4     19888     Rp1    chr1
5     20671   Sox17    chr1
27174 more rows ...
```


## 2. 数据预处理

#### 2.1 转换count数据为CPM值
差异表达分析之前，由于不同的测序深度会导致counts数目不同，所以为了消除这种由建库引起的基因之间的差异，我们需要将数据归一化。归一化常见的方法有： `CPM`, `log-CPM`, `RPKM` 和 `FPKM`  。本例中使用CPM值（这已足够）。
+ 本例中使用 **`edgeR`** 包中的 **`cpm`** 函数来完成数据的归一化，log转换的时候取0.25作为先验来避免取到0 . （**edgeR** 包也提供了 `rpkm` 函数 来计算 RPKM值）
```r
## convert counts to CPM and log-CPM
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE) # 取log值
```
#### 2.2 过滤表达量太低的基因

```r
## Removing genes that are lowly expressed
# 首先查看那些在所有样本中均为0的基因数目
table(rowSums(x$counts==0)==9) # 可以发现有5153个基因在所有样本中表达量均为0
```
+ 在任何情况下，生物学水平上没有表达的基因都应该被丢弃，这样可以将基因的子集缩小到那些感兴趣的基因，并在做差异表达时减少下游分析所用于检验的基因数量。通过检验log-CPM值可以看出，每个样本中都有很大比例的基因未表达或低表达(如图a所示)。
+ 以CPM为1(相当于log-CPM值为0)作为分界线阈值来检查基因是低表达或者高表达，如果其表达高于此阈值，则视为表达，否则为低表达。基因必须在至少一组中表达(或在整个实验中至少三个样本中表达)，以备后续分析。
+ 在这里，CPM值为1意味着如果一个基因在测序深度最低的样本中至少有20个计数(JMS9-P8c，文库大小约为2600万)。在测序深度最大的样本中至少有76个计数(JMS8-3，库大小约为。7600万)。如果测得的reads是外显子的而不是整个基因的，或实验的测序深度较低，则可以考虑较低的CPM阈值。

```r
# 取出至少在三个样本中cpm值均大于1的基因
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
```

