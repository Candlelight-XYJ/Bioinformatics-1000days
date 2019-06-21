setwd("E:/GitHub/Bioinformatics-1000days/day16~day22-RNA-seq123 with limma edgeR and Glimma/")
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

## Use readDEG read counts 将9个样本转为DGE对象
x <- readDGE(file.path(".", files), columns=c(1,3))
class(x)
dim(x)

## Organising sample information 设定sample的名字
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames

## 设置sample信息
colnames(x) <- samplenames
# 分组
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane # 添加lane的信息
x$samples # 此时sample信息中就会多出分组和细胞类型两列

## using the Mus.musculus package 
## to retrieve associated gene symbols and chromosome information

geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID") # 提取symbol和chr, 以entrezid作为map id的源头
head(genes)
## !duplicated replicate genes
genes <- genes[!duplicated(genes$ENTREZID),]
## add genes 将提取好的基因信息添加入我们的DGEList中genes这个二级数据框
x$genes <- genes
x

## convert counts to CPM and log-CPM
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE) # 取log值

## Removing genes that are lowly expressed
# 首先查看那些在所有样本中均为0的基因数目
table(rowSums(x$counts==0)==9) # 可以发现有5153个基因在所有样本中表达量均为0
# 取出至少在三个样本中cpm值均大于1的基因
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE] ## 代码啥意思
dim(x)

## 绘图展示基因过滤前后log-cpm值的分布情况
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired") ## 配置绘图调色盘的主题，paired 是 qualitative palettes 中的一个颜色配置
# col <- brewer.pal(nsamples, "Pastel1") # 尝试使用其它颜色
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="") # 首先对第一列的第一个样本绘图
title(main="A. Raw data", xlab="Log-cpm") # 加入title
abline(v=0, lty=3) # 在0坐标处添加分隔虚线
for (i in 2:nsamples){ # 批量在同一画布上绘出其它样本的cpm分布情况
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n") # 添加图例
## 开始对过滤后的lcpm绘图,函数功能注释基本同上
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="") 
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

## 使用TMM标准化基因表达分布
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

## 可视化准备
# 可视化数据预处理
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5
# 绘制箱线图
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
#> [1] 0.05472223 6.13059440 1.22927355 1.17051887 1.21487709 1.05622968
#> [7] 1.14587663 1.26129350 1.11702264
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

## 非监督聚类
# 展示样本之间的差异性与相似性
# 同一个处理的多个重复，没啥实验误差的话一般会聚在一起
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2)) # set huabu
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")
# Glimma 提供的绘MDS图函数
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)


## 差异表达分析
# First,建立分组信息
design <- model.matrix(~0+group+lane) # 设置分组矩阵
colnames(design) <- gsub("group", "", colnames(design)) # 去掉列名中的group
design

# Second,use makeContrasts function 建立比较信息
contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP, # Basal 和 LP 比较
  BasalvsML = Basal - ML, # Basal 和 ML 比较
  LPvsML = LP - ML, # LP 和 ML 比较
  levels = colnames(design)) 
contr.matrix

# Third, Removing heteroscedascity from count data
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE) # Transform RNA-Seq Data Ready for Linear Modelling
v

vfit <- lmFit(v, design) # Fit linear model for each gene given a series of arrays
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) # Given a linear model fit to microarray data, compute estimated coefficients and standard errors for a given set of contrasts.
efit <- eBayes(vfit) # Empirical Bayes Statistics for Differential Expression
plotSA(efit, main="Final model: Mean-variance trend")

# result
summary(decideTests(efit))

# treat 类似于 ebayes
# When the number of DE genes is large, 
# treat is often useful for giving preference to 
# larger fold-changes and for prioritizing genes 
# that are biologically important
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

# Fourth, plot Venn 
# 取出BasalvsLP和BasalvsML 这两组比较中的共同差异基因
de.common <- which(dt[,1]!=0 & dt[,2]!=0) 
# 查看共同差异基因数目
length(de.common)
# 查看前20个基因symbol
head(tfit$genes$SYMBOL[de.common], n=20)
# 绘制Venn图
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

# Fifth , output result
write.fit(tfit, dt, file="results.txt")

## Examining individual DE genes from top to bottom

# 使用 topTreat() 将差异基因按padj,logFC,log-CPM,t值 从小到大排序
# n=Inf 表示选取所有基因参与排序
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf) # coef在此处代表选取的比对的组别
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)
head(basal.vs.ml)

## 差异基因结果可视化
# plot MD
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))
# plot MD using GLimma 这个包绘制出来的是动态的MD图
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=x$counts, groups=group, launch=FALSE)


# plot heatmap
library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

