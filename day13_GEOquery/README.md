## 学习章节
> https://bioconductor.github.io/BiocWorkshops/public-data-resources-and-bioconductor.html


#### 学习目标
+ 学会使用Bioconductor包去获取并操作公共数据库中的数据
+ 包括GEO（Gene Expression Omnibus ），SRA（Sequence Read Archive）,GDC（Genomic Data Commons ），存储在Bioconductor上的宏基因组数据，药物基因组数据（ PharmacoDB）以及癌症基因组数据

#### 需要预先准备的R包
```r
library(GEOquery)
library(GenomicDataCommons)
library(SRAdbV2)
library(curatedTCGAdata)
library(curatedMetagenomicData)
library(HMP16SData)
library(PharmacoGx)
``` 

---

## 1. GEOquery
NCBI Gene Expression Omnibus (GEO)是一个公共存储库，它存储了大量的高通量的实验数据。这些数据包括基于芯片的单通道和双通道实验数据，例如：mRNA表达数据、基因组DNA数据和蛋白质丰度数据。它也包含了非芯片技术产生的数据，例如基因表达序列分析(SAGE)数据、质谱蛋白质组数据和高通量测序数据。 **`GEOquery`** 是Bioconductor中用于获取GEO数据的R包。

#### 1.1 Overview of GEO
GEO数据库由4个板块构成，前三个（Sample, Platform, and Series）由用户提供；后一个（datasets）数据集，由GEO工作人员根据用户提供的数据进行整合管理。

注：GDS已经停止提供了

#### 1.2 GEOquery的使用案例：MDS plot of cancer data
使用GEOquery包中的 **`getGEO`**函数即可快速获取GEO的数据
+ **本次案例数据来源**：https://doi.org/10.1158/1055-9965.EPI-17-0461
> 背景:肿瘤微环境是影响肿瘤免疫治疗反应的重要因素。为了进一步了解肿瘤如何影响局部免疫系统，我们分析了匹配正常组织和肿瘤组织之间的免疫基因表达差异。方法:我们分析了来自实体癌症和分离免疫细胞群的公开的和新的基因表达数据。我们还确定了CD8、FoxP3免疫组化和我们的基因签名之间的相关性。结果:调节T细胞(Tregs)是正常组织和肿瘤组织免疫基因表达差异的主要驱动因素之一
+ **本次案例涉及练习**
  + 使用GEOquery获取公共组学数据
  + 将公共组学数据转换为 **`SummarizedExperiment`** 对象
  + 可视化这些公共数据

```r
# download data from GEO
gse = getGEO("GSE103512")[[1]]

# convert the old ExpressionSet structure to the newer SummarizedExperiment
library(SummarizedExperiment)
se = as(gse, "SummarizedExperiment")

# Examine two variables of interest, cancer type and tumor/normal status.
with(colData(se),table(`cancer.type.ch1`,`normal.ch1`))

# Filter gene expression by variance to find most informative genes.
sds = apply(assay(se, 'exprs'),1,sd)
dat = assay(se, 'exprs')[order(sds,decreasing = TRUE)[1:500],]

# Perform multidimensional scaling and prepare for plotting
# make a data.frame before plotting
mdsvals = cmdscale(dist(t(dat)))
mdsvals = as.data.frame(mdsvals)
mdsvals$Type=factor(colData(se)[,'cancer.type.ch1'])
mdsvals$Normal = factor(colData(se)[,'normal.ch1'])
head(mdsvals)

# do the plot
library(ggplot2)
ggplot(mdsvals, aes(x=V1,y=V2,shape=Normal,color=Type)) + 
  geom_point( alpha=0.6) + theme(text=element_text(size = 18))
```

#### 1.3 Accessing Raw Data from GEO
> 背景：NCBI GEO accepts (but has not always required) raw data such as .CEL files, .CDF files, images, etc. It is also not uncommon for some RNA-seq or other sequencing datasets to supply only raw data (with accompanying sample information, of course), necessitating Sometimes, it is useful to get quick access to such data

我们可以使用 **`getGEOSuppFiles（）`** 函数获取raw data（以GEO的Accession number作为参数，例如:GSE12387等等）。 默认情况下，这个函数会在当前的工作目录下自动创建一个文件夹，来存储用户选择下载的raw data。









