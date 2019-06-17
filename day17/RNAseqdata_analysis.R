setwd("E:/GitHub/Bioinformatics-1000days/day17/")
## install RNASeq123
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



suppressPackageStartupMessages({
  library(limma)
  library(Glimma)
  library(edgeR)
  library(Mus.musculus)
})

## read data
read.delim(file.path(".", files[1]), nrow=5)

## Use readDEG read counts
x <- readDGE(file.path(".", files), columns=c(1,3))
class(x)
dim(x)

## Organising sample information

