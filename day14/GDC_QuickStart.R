## report bugs
bug.report(package='GenomicDataCommons') 

## installation
install.packages('BiocManager')
BiocManager::install('GenomicDataCommons')
library(GenomicDataCommons)

## check this connectivity and functionality
GenomicDataCommons::status()


## filtering finds gene expression files quantified as raw counts using HTSeq from ovarian cancer patients
ge_manifest = files() %>%
  filter( ~ cases.project.project_id == 'TCGA-OV' &
            type == 'gene_expression' &
            analysis.workflow_type == 'HTSeq - Counts') %>%
  manifest()

## download data
fnames = lapply(ge_manifest$id[1:2],gdcdata)

## Metadata queries
expands = c("diagnoses","annotations",
            "demographic","exposures")
projResults = projects() %>%
  results(size=10)
str(projResults,list.len=5)
names(projResults)
