## 学习章节
https://bioconductor.github.io/BiocWorkshops/r-and-bioconductor-for-everyone-an-introduction.html#working-with-summarized-experimental-data

[toc]

## 1. Working with summarized experimental data
#### 1.1 简介
本章主要学习SummatizedExperiment包和`SummarizedExperiment`对象
`SummarizedExperiment`对象具有类似于矩阵的性质，我们可以通过行和列，对它取子集。
来自于SummarizedExperiment对象实验的数据assay()，它的行代表我们感兴趣的特征(例如基因)，列代表每个样本，（矩阵中的每个值可能代表每个基因的在不同样本中的表达量)
![SummarizedExperiment](https://github.com/Candlelight-XYJ/Bioinformatics-1000days/blob/master/day2/SummarizedExperiment.png)

#### 1.2 构建SummarizedExperiment对象
> **数据介绍**
> 包含有8个样本，数据由RNA-seq实验产生，主要是用于观察4个人的平滑肌细胞系对地塞米松治疗的情况
> 我们可以使用函数` browseVignettes("airway")`查看关于这个数据集和实验的详细描述
```r
## input data
fname <- file.choose() # airway_colData.csv
fname
## set the first column of the data to be treated as row names(将第一列作为数据的row-names)
colData <- read.csv(fname, row.names = 1)
colData
```
这组数据来源于Short Read Archive，包含`SampleName`,`Run`,`Experiment`,`Sampel`,`BioSample`这些列，另外我们还需要添加以下的列：
+ Cell:所使用的细胞系，本数据使用了4个细胞系
+ dex:这个样本是否添加了地塞米松
+ albut:二次治疗，我们可以忽略
+ avgLength:本次实验中，每个样本的RNA-seq的reads的平均长度

#### 1.3 Assay data
现在导入assay数据
```r
## importing the assay data from the file “airway_counts.csv”
fname <- file.choose() # airway_counts.csv
fname

counts <- read.csv(fname, row.names=1)
## coerce data.frame() to matrix using as.matrix()
counts <- as.matrix(counts)
## We see the dimensions and first few rows of the counts matrix
dim(counts)
#> [1] 33469 8
head(counts)
```
数据解释
+ 以基因ENSG00000000003为例，样本SRR1039508 有679 个reads，覆盖了它；样本SRR1039509 有448个reads覆盖了它。

#### 1.4 Creating a SummarizedExperiment object
```r
## Attach the SummarizedExperiment library to our R session
library("SummarizedExperiment")
## Use the SummarizedExperiment() function to coordinate the assay and column data
## 校准数据
se <- SummarizedExperiment(assay = counts, colData = colData)
se
## use subset() on SummarizedExperiment to create subsets of the data in a coordinated way
## 取出数据中的子集，注意由于SummarizedExperiment是个二维矩阵，所以我们对他的操作也是基于二维的
subset(se, , dex == "trt")
## use assay() to extract the count matrix, 
## colSums() to calculate the library size (total number of reads overlapping genes in each sample)
## colSums()计算每个样本中覆盖了所有基因的reads总数
colSums(assay(se))
## 
se$lib.size <- colSums(assay(se))
colData(se)
```
## 2. 下游分析 Down-stream analysis
使用R包`DESeq2`来进行下游分析
```r
## Down-stream analysis
library("DESeq2")
## including cell line as a covariate, 
## and dexamethazone treatment as the main factor that we are interested in
## 构建dds数据集
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds
## performs advanced statistical analysis on the data in the dds object
## 进行统计分析
dds <- DESeq(dds)
## A table summarizing measures of differential expression can be extracted from the object
## 使用results查看差异分析结果
results(dds)
```