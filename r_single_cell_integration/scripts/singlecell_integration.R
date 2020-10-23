renv::restore()
BiocManager::install('Seurat')

library(tidyverse)
library(batchelor)
library(scran)
library(Seurat)
library(scater)
library(DropletUtils)
library(patchwork)

sce <- read10xCounts(
    c(v2 = 'data/pbmc_1k_v2/',
      v3 = 'data/pbmc_1k_v3/'), col.names = TRUE)




# Exploring single cell experiment object (sce)

sce
dim(sce)
str(sce)
slotNames(sce)
metadata(sce) # Miscellaneous list of extra metadata

colnames(colData(sce)) # Cell metadata
colData(sce)

colnames(rowData(sce)) # Gene metadata
rowData(sce)

assayNames(sce)
assay(sce)

# QC=======

is.mito <- grepl("^MT-", rowData(sce)$Symbol)
table(is.mito)

qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
#quickperCellQC uses the isoutlier function
#convenience function to determine which values in a numeric vector are outliers based on the median absolute deviation (MAD)
table(filtered)
sce <- sce[, !filtered$discard]

# Normalisation pre-filtering=====

sce <- logNormCounts(sce)
range(assay(sce,"logcounts")) # 0.00000 12.28542

decomposed_var <- modelGeneVar(sce)

HVG <- getTopHVGs(decomposed_var, prop=0.1)



# Quick exploration=======
# These datasets are pre-filtered

set.seed(123)
sce <- runPCA(sce,
              name='RNA_PCA',
              ntop=Inf, # Use all the HVGs I give you
              subset_row=HVG)

str(reducedDims(sce)$RNA_PCA)
colData(sce)

plotReducedDim(sce,
               dimred = "RNA_PCA",
               ncomponents = 2,
               colour_by = "Sample")

#Can we decide whether a linear batch correction is sufficient just off PCA?
#PRObaby not

# Scree plot to choose PCs for UMAP


percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, xlab="PC", ylab="Variance explained (%)")


sce <- runUMAP(sce,
               name="UMAP",
               dimred="RNA_PCA",
               n_dimred=20,
               colour_by = "Sample"
)

#colour_by = "Sample"
# 20 Looks good on the scree plot


plotReducedDim(sce,
               dimred = "UMAP",
               ncomponents = 2)

#Integrate v2 and v3 datasets =======

#We use the rescaleBatches() function from the batchelor package to remove the batch effect. This is roughly equivalent to applying a linear regression to the log-expression values per gene, with some adjustments to improve performance and efficiency. For each gene, the mean expression in each batch is scaled down until it is equal to the lowest mean across all batches

sce_v2 <- sce[,sce$Sample == 'v2']
dim(sce_v2)
sce_v3 <- sce[,sce$Sample == 'v3']
dim(sce_v3)


rescaled <- rescaleBatches(sce_v2, sce_v3)
rescaled
dim(rescaled)
colData(rescaled)

set.seed(123) # To ensure reproducibility of IRLBA:NOTE: using 123
rescaled <- runPCA(rescaled,
                   subset_row=HVG,
                   exprs_values="corrected")

reducedDim(rescaled)
plotReducedDim(rescaled, dimred = 'PCA',
               ncomponents = 2,
               colour_by = "cluster")

#in this case, the linear batch correction method employed in rescaledBatches
#wasn't right to fully correct for the differences across all 'clusters'
#eg the right most cluster
#not right method to correct all cells in correct linear regression
#However these is no linear method to FULLY correct for all differences among batches, so this isn't a terrible outcome

#One way of checking batch correct is to use SNN
#If integration works well, then cells from all batches should be present in
#all clusters
snn.gr <- buildSNNGraph(rescaled, use.dimred="PCA")
clusters.resc <- igraph::cluster_walktrap(snn.gr)$membership
colData(rescaled)$cluster_walktrap_SNN <- as.factor(clusters.resc)
tab.resc <- table(Cluster=clusters.resc, Batch=rescaled$batch)
tab.resc

plotReducedDim(rescaled, dimred = 'PCA',
               ncomponents = 2,
               colour_by = "cluster_walktrap_SNN",
               shape_by = 'batch')

rescaled <- runUMAP(rescaled,
               name="UMAP",
               dimred="PCA",
               n_dimred=20)

p1 <- plotReducedDim(rescaled, dimred = 'UMAP',
               ncomponents = 2,
               colour_by = "cluster_walktrap_SNN",
               shape_by = 'batch')


p2<- plotReducedDim(rescaled, dimred = 'UMAP',
               ncomponents = 2,
               colour_by = "batch")

#plot UMAPs together using patchwork
p1 + p2

#haven't checked scree plot, assume 20 is fine from previous check

#how many cells from each batch are in each cluster


#Performing MNN correction=======
