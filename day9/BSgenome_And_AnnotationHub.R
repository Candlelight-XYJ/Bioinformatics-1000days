##############
## BSgenome ##
##############

# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# get a list of available BSgenome packages
library(BSgenome)
head(available.genomes())

#  load and inspect a BSgenome package
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens

# get data by sequence ,or by passing in a GRanges object, to get just a region
getSeq(Hsapiens, "chr1")

getSeq(Hsapiens, gns["5467",])

# Get the sequences for all transcripts of the TP53 gene
browseVignettes("BSgenome")


####################
##  AnnotationHub ##
####################

## simple use
library(AnnotationHub)
hub <- AnnotationHub()
#> snapshotDate(): 2018-06-27
hub

names(mcols(hub))


## AnnotationHub Data providers
unique(hub$dataprovider)

## AnnotationHub Data classes
unique(hub$rdataclass)

## AnnotationHub Species
head(unique(hub$species))
length(unique(hub$species))

## Data sources
unique(hub$sourcetype)

## AnnotationHub query
qry <- query(hub, c("granges","homo sapiens","ensembl"))
qry
qry$sourceurl

## Selecting AnnotationHub resource
whatIwant <- qry[["AH50377"]]

GRCh38TxDb <- makeTxDbFromGRanges(whatIwant)
GRCh38TxDb



