#############################
## step1 download datasets ##
#############################

library(GEOquery)
eSet <- getGEO("GSE54238",destdir=".",AnnotGPL = F, getGPL = F)

## select exp
expressMatrix <- exprs(eSet[[1]])
## select pheno data
pheno <- pData(eSet[[1]])
## group info
library(stringr)
group <- str_split()


