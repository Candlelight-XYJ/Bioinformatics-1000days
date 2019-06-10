## 学习章节
https://bioconductor.github.io/BiocWorkshops/public-data-resources-and-bioconductor.html#genomicdatacommons
[toc]

## 1. GenomicDataCommons
+ 关于GDC介绍
> The National Cancer Institute’s (NCI’s) Genomic Data Commons (GDC) is a data sharing platform that promotes precision medicine in oncology. It is not just a database or a tool; it is an expandable knowledge network supporting the import and standardization of genomic and clinical data from cancer research programs. The GDC contains NCI-generated data from some of the largest and most comprehensive cancer genomic datasets, including The Cancer Genome Atlas (TCGA) and Therapeutically Applicable Research to Generate Effective Therapies (TARGET). For the first time, these datasets have been harmonized using a common set of bioinformatics pipelines, so that the data can be directly compared. As a growing knowledge system for cancer, the GDC also enables researchers to submit data, and harmonizes these data for import into the GDC. As more researchers add clinical and genomic data to the GDC, it will become an even more powerful tool for making discoveries about the molecular basis of cancer that may lead to better care for patients

GDC数据库的数据结构较为复杂(如下图所示)，详细可访问:https://gdc.cancer.gov/developers/gdc-data-model/gdc-data-model-components

[gdc](https://github.com/Candlelight-XYJ/Bioinformatics-1000days/blob/master/day14/gdc.png)

数据模型被编码为一个属性图。这个图中的节点实体代表 **`项目、病例、诊断、档案(各种类型)以及注释`** . 实体之间的关系由边相连。 节点和边都可能提供实例较为详细信息的属性。
＞　The GDC API exposes these nodes and edges in a somewhat simplified set of RESTful endpoints

## 2. Quick start
反馈R包的bug，可以通过： **`bug.report(package='GenomicDataCommons')`** 的方式来转入GitHub的Issue页进行提交

#### 2.1 安装包&检查数据库连接状态
+ 首先安装GDC包
```r
## report bugs
bug.report(package='GenomicDataCommons') 

## installation
install.packages('BiocManager')
BiocManager::install('GenomicDataCommons')
library(GenomicDataCommons)
```
+ 然后检查GDC包和数据库的网络连接状态
要保证NCI的GDC的API接口可用才行
可以使用 **`status()`** 函数来检查包的网络连接情况和功能
```r
$commit
[1] "3e22a4257d5079ae9f7193950b51ed9dfc561ed1"

$data_release
[1] "Data Release 17.0 - June 05, 2019"

$status
[1] "OK"

$tag
[1] "1.21.0"

$version
[1] 1
```

#### 2.2 查找数据
+ 一个通过R查找GDC原始数据的典型示例（下载来自于卵巢癌患者，测序平台为 `HTSeq` 经过定量后的raw counts基因表达数据）
```r
## filtering finds gene expression files quantified as raw counts using HTSeq from ovarian cancer patients
ge_manifest = files() %>%
    filter( ~ cases.project.project_id == 'TCGA-OV' &
                type == 'gene_expression' &
                analysis.workflow_type == 'HTSeq - Counts') %>%
    manifest()
```

#### 2.3 下载数据
上一节代码查找到了379个满足条件的基因表达数据，这一节来看如何下载这些被查找到的数据。
+ 在许多情况下，使用多个进程进行下载可以显著加快传输速度
+ 第一次下载数据的时候，R会要求建立一个数据下载缓存的目录（通过 **`?gdc_cache`** 查看细节），下载的结果将存储于缓存目录
+ 之后对相同文件的访问将直接从缓存目录中进行，从而避免了多次下载同一个数据

```r
fnames = lapply(ge_manifest$id[1:2],gdcdata)
## Would you like to create a GDC Cache directory at C:\Users\ADMINI~1\AppData\Local\GenomicDataCommons\GenomicDataCommons\Cache 

## 1: Yes
## 2: No
## Selection: 1
  |============================================================================================| 100%

  |============================================================================================| 100%
```

#### 2.4 Metadata queries
GenomicDataCommons可以获取来自NCI GDC的临床，人口，生物特异性以及注释的信息

```r
expands = c("diagnoses","annotations",
             "demographic","exposures")
projResults = projects() %>%
    results(size=10)
str(projResults,list.len=5)
# List of 9
# $ dbgap_accession_number: chr [1:10] "phs001287" "phs001374" "phs001628" "phs000466" ...
# $ disease_type          :List of 10
#  ..$ CPTAC-3              : chr "Adenomas and Adenocarcinomas"
#  ..$ VAREPOP-APOLLO       : chr [1:2] "Epithelial Neoplasms, NOS" "Squamous Cell Neoplasms"
#  ..$ BEATAML1.0-CRENOLANIB: chr "Myeloid Leukemias"
#  ..$ TARGET-CCSK          : chr "Clear Cell Sarcoma of the Kidney"
#  ..$ TARGET-NBL           : chr "Neuroblastoma"
#  .. [list output truncated]
# $ releasable            : logi [1:10] FALSE FALSE FALSE TRUE TRUE FALSE ...
# $ released              : logi [1:10] TRUE TRUE TRUE TRUE TRUE TRUE ...
# $ state                 : chr [1:10] "open" "open" "open" "open" ...
#  [list output truncated]
# - attr(*, "row.names")= int [1:10] 1 2 3 4 5 6 7 8 9 10
# - attr(*, "class")= chr [1:3] "GDCprojectsResults" "GDCResults" "list"

names(projResults)
#[1] "dbgap_accession_number" "disease_type"           "releasable"            
#[4] "released"               "state"                  "primary_site"          
#[7] "project_id"             "id"                     "name"   


```






