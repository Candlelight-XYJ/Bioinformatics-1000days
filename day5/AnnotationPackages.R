library(TxDb.Hsapiens.UCSC.hg38.knownGene)

## the genomic position of all genes
gns <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gns


## the genomic position of all transcripts by gene
txs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene)
txs


## use ` [ ` function to subset GRanges and GRangesLists Objects
txs[txs %over% gns[1:2,]]