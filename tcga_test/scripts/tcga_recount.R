#Downloading TCGA data with recount package 
#RangeSummarizedExperiment object 


#Install==== 

## Install packages from Bioconductor
install.packages("BiocManager")
BiocManager::install(c(
  "recount", "GenomicRanges", "DESeq2", "ideal", "regionReport",
  "clusterProfiler", "org.Hs.eg.db", "gplots", "derfinder",
  "rtracklayer", "GenomicFeatures", "bumphunter", "derfinderPlot",
  "devtools"))
```

#Library==== 
library(recount)
library(GenomicRanges)
library(DESeq2)
library(ideal)
library(regionReport)
library(clusterProfiler)
library(org.Hs.eg.db) #failed to install 
library(ggplot2)
library(derfinder)
library(rtracklayer)
library(GenomicFeatures)
library(bumphunter)
library(derfinderPlot)



#recount package=== 
#Find project of interest: 

tcga_paad_info <- abstract_search('TCGA-PAAD')
