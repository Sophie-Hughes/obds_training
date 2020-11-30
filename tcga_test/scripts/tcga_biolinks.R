#TCGABiolinks to download TCGA data# 

#Install===  Remove this from script - or have in seperate script (r_install)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

BiocManager::install(c('ggplot2', 'recount'))


#Library=== 
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(DESeq2)



#Define a list of samples to query and download:
  ## these are barcodes are taken from TCGA-PAAD, adenocarcinoma only 

listSamples <- c('TCGA-HZ-A9TJ', 'TCGA-HZ-7289', 'TCGA-US-A77E','TCGA-YH-A8SY',
                 'TCGA-HZ-A77Q', 'TCGA-RB-A7B8', 'TCGA-HZ-A4BK', 'TCGA-HZ-A77P',
                 'TCGA-3A-A9IV', 'TCGA-2J-AABH', 'TCGA-F2-6879', 'TCGA-F2-A44G',
                 'TCGA-3A-A9IS', 'TCGA-H6-8124', 'TCGA-2L-AAQM', 'TCGA-2J-AABI',
                 'TCGA-2J-AAB4', 'TCGA-3A-A9IN', 'TCGA-H6-A45N', 'TCGA-H8-A6C1')

#GDC Query====
# Query platform Illumina HiSeq with a list of barcode 
paad_query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples, 
                  legacy = TRUE)

#GDC Download====
# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(paad_query)

#GDC Prepare====
# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
PDACRnaseqSE <- GDCprepare(paad_query)

PDACMatrix <- assay(PDACRnaseqSE,"raw_count") # or PDACMatrix <- assay(PDACRnaseqSE,"raw_count")

# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run:
PDACRnaseq_CorOutliers <- TCGAanalyze_Preprocessing(PDACRnaseqSE)

rowRanges(PDACRnaseqSE)
colData(PDACRnaseqSE)

