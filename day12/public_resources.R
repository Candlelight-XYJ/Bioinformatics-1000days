BiocManager::install(c("GenomicDataCommons","SRAdbV2","curatedTCGAdata",
                       "curatedMetagenomicData","HMP16SData","PharmacoGx"))

library(GEOquery)
library(GenomicDataCommons)
library(SRAdbV2)
library(curatedTCGAdata)
library(curatedMetagenomicData)
library(HMP16SData)
library(PharmacoGx)

######################################
## case1 : MDS plot of cancer data  ##
######################################

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



