######################################
## step1 download datasets GSE54238 ##
######################################

library(GEOquery)
eSet <- getGEO("GSE54238",destdir=".",AnnotGPL = F, getGPL = F)

## select exp
expressMatrix <- exprs(eSet[[1]])
## select pheno data
pheno <- pData(eSet[[1]])
## group info
library(stringr)
group <- paste0(str_split(pheno$source_name_ch1,' ',n=2,simplify = T))
## save data
save(group,expressMatrix,file="GSE54238.Rdata")

######################################
## step1 download datasets GSE14520 ##
######################################
eSet2 <- getGEO("GSE14520",destdir=".",AnnotGPL = F, getGPL = F)
save(eSet2,file="GSE14520.Rdata")
load("GSE14520.Rdata")
## select exp
expressMatrix2 <- exprs(eSet[[1]])
## select pheno data
pheno2 <- pData(eSet[[1]])
## group info
library(stringr)
group2 <- paste0(str_split(pheno$source_name_ch1,' ',n=2,simplify = T))
## save data
save(group2,expressMatrix2,file="GSE14520.Rdata")

######################################
## step2  ##
######################################





