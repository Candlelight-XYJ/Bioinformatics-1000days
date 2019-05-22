## input data
fname <- file.choose()  # airway_colData.csv
fname
## set  the first column of the data to be treated as row names
colData <- read.csv(fname, row.names = 1)
colData
## importing the assay data from the file “airway_counts.csv”
fname <- file.choose()  # airway_counts.csv
fname

counts <- read.csv(fname, row.names=1)
## coerce data.frame() to matrix using as.matrix()
counts <- as.matrix(counts)
## We see the dimensions and first few rows of the counts matrix
dim(counts)
#> [1] 33469     8
head(counts)


## Attach the SummarizedExperiment library to our R session
library("SummarizedExperiment")
## Use the SummarizedExperiment() function to coordinate the assay and column data
se <- SummarizedExperiment(assay = counts, colData = colData)
se

## use subset() on SummarizedExperiment to create subsets of the data in a coordinated way
subset(se,, dex == "trt")

## use assay() to extract the count matrix, 
## colSums() to calculate the library size (total number of reads overlapping genes in each sample)
## colSums()计算每个样本中覆盖了所有基因的reads总数
colSums(assay(se))


## 
se$lib.size <- colSums(assay(se))
colData(se)


## Down-stream analysis
library("DESeq2")
## including cell line as a covariate, 
## 构建dds数据集
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds
## performs advanced statistical analysis on the data in the dds object
## 进行统计分析
dds <- DESeq(dds)
## A table summarizing measures of differential expression can be extracted from the object
## 使用results查看差异分析结果
results(dds)

## 接下来再使用可视化等方法对此结果处理即可

