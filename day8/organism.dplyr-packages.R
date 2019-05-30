#################################
## organism packages exercises ##
#################################


## Get all the GO terms for BRCA1


## What gene does the UCSC transcript ID uc002fai.3 map to?
  
## How many other transcripts does that gene have?
  
## Get all the transcripts from the hg19 genome build, along with their Ensembl gene ID, UCSC transcript ID and gene symbol


####################
## Organism.dplyr ##
####################
library(Organism.dplyr)
src <- src_organism(dbpath = hg38light())
src

## Get promoters from a TxDb object
options(width = 120)
promoters(src)

## Extract a table from the underlying database
tbl(src, "id")


## Make a complex query between tables in the underlying database
inner_join(tbl(src, "id"), tbl(src, "ranges_gene")) %>%
  filter(symbol %in% c("ADA", "NAT2")) %>%
  dplyr::select(gene_chrom, gene_start, gene_end,
                gene_strand, symbol, alias, map)

##############################
## Organism.dplyr exercises ##
##############################
## How many supported organisms are implemented in Organism.dplyr?
browseVignettes("Organism.dplyr")
supportedOrganisms()

## Display the ensembl Id and genename description for symbol “NAT2”
# look at all avaliable tables
src_tbls(src)
## get ensembl Id and genename description
  tbl(src, "id") %>%
  filter(symbol == "NAT2") %>%
  dplyr::select( ensembl,genename)

## Show all the alias for “NAT2” in the database
tbl(src, "id") %>%
  filter(symbol == "NAT2") %>%
  dplyr::select(alias)

## Display Gene ontology (GO) information for gene symbol “NAT2”
inner_join(tbl(src, "id"), tbl(src, "id_go")) %>%
  filter(symbol == "NAT2") %>%
  dplyr::select(entrez, ensembl, symbol, go, evidence, ontology)




