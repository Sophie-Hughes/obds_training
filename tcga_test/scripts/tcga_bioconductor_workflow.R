# TCGA workflow taken from Noushmehr et al, 2020 
# https://www.bioconductor.org/packages/release/workflows/vignettes/TCGAWorkflow/inst/doc/TCGAWorkflow.html 


#Installation======  
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

#Library Installation===
library(TCGAbiolinksGUI.data)
library(dplyr)
library(DT)

#Searching GDC data for download====== 

getGDCdisease <- function(){
  projects <- TCGAbiolinks:::getGDCprojects()
  projects <- projects[projects$id != "FM-AD",]
  disease <-  projects$project_id
  idx <- grep("disease_type",colnames(projects))
  names(disease) <-  paste0(projects[[idx]], " (",disease,")")
  disease <- disease[sort(names(disease))]
  return(disease)
}

data(GDCdisease)
DT::datatable(as.data.frame(GDCdisease))
