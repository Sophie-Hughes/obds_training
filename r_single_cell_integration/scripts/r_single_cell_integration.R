renv::restore()
BiocManager::install("batchelor")
BiocManager::install("Seurat")
library(tidyverse)
library(batchelor)
library(scran)
library(Seurat)
library(scater)
library(DropletUtils)

# https://osca.bioconductor.org/integrating-datasets.html#linear-regression
# Load in data/create sce -----

#Loading data and creating single cell experiment object
sce <- read10xCounts(
    c(v2="data/pbmc_1k_v2/",
      v3='data/pbmc_1k_v3/'),
    col.names = TRUE)


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

# QC -----
is.mito <- grepl("^MT-", rowData(sce)$Symbol)
table(is.mito)
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
# quickPerCellQC uses the isOutlier function
# # Convenience function to determine which values in a numeric vector are
# # outliers based on the median absolute deviation (MAD).

filtered
table(filtered)
sce <- sce[, !filtered$discard]


# Normalisation post-filtering ----

sce <- logNormCounts(sce)
range(assay(sce,"logcounts")) # 0.00000 11.28685

decomposed_var <- modelGeneVar(sce)
HVG <- getTopHVGs(decomposed_var, prop=0.1)

# Quick exploration -----
# These datasets are pre-filtered

set.seed(123)
sce <- runPCA(sce,
              name='PCA',
              ntop=Inf, # Use all the HVGs I give you
              subset_row=HVG)
str(reducedDims(sce)$RNA_PCA)
colData(sce)
plotReducedDim(sce,
               dimred = "PCA",
               ncomponents = 2,
               colour_by = "Sample")
# Can we decide whether a linear batch correction is sufficient just off PCA?
# Probably not

# Scree plot to choose PCs for UMAP

percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, xlab="PC", ylab="Variance explained (%)")

sce <- runUMAP(sce,
               name="UMAP",
               dimred="PCA",
               n_dimred=20)
# 20 Looks good on the scree plot

plotReducedDim(sce,
               dimred = "UMAP",
               ncomponents = 2)

# Integrate v2 and v3 datasets -----
# We use the rescaleBatches() function from the batchelor package to remove
# the batch effect. This is roughly equivalent to applying a linear regression
# to the log-expression values per gene, with some adjustments to improve
# performance and efficiency. For each gene, the mean expression in each batch
# is scaled down until it is equal to the lowest mean across all batches.

sce_v2 <- sce[,sce$Sample=="v2"]
dim(sce_v2)
sce_v3 <- sce[,sce$Sample=="v3"]
dim(sce_v3)

rescaled <- rescaleBatches(sce_v2, sce_v3)
rescaled
dim(rescaled)
colData(rescaled)

set.seed(123) # To ensure reproducibility of IRLBA. NOTE: USING 123
rescaled <- runPCA(rescaled,
                   subset_row=HVG,
                   exprs_values="corrected")

str(reducedDim(rescaled))
plotReducedDim(rescaled,
               dimred = "PCA",
               ncomponents = 2,
               colour_by = "batch")
# In this case, the linear batch correction method employed in rescaleBatches
# was insufficient to fully correct for the differences across all 'clusters'
# (eg the right-most cluster)
# However, there is no linear method to FULLY correct for all differences among
# batches, so this isn't a terrible outcome.

# One way of checking batch correction is to use SNN
# If integration works well, then cells from all batches should be present in
# all clusters
snn.gr <- buildSNNGraph(rescaled, use.dimred="PCA")
clusters.resc <- igraph::cluster_walktrap(snn.gr)$membership
colData(rescaled)$cluster_walktrap_SNN <- as.factor(clusters.resc)
tab.resc <- table(Cluster=clusters.resc, Batch=rescaled$batch)
tab.resc

plotReducedDim(rescaled,
               dimred = "PCA",
               ncomponents = 2,
               colour_by = "cluster_walktrap_SNN",
               shape_by = "batch")

rescaled <- runUMAP(rescaled,
               name="UMAP",
               dimred="PCA",
               n_dimred=20)
# Haven't checked scree plot, assume 20 is fine from previous check
library(patchwork)

colData(rescaled)$batch <- as.factor(colData(rescaled)$batch)

p1 <- plotReducedDim(rescaled,
               dimred = "UMAP",
               ncomponents = 2,
               colour_by = "cluster_walktrap_SNN")
rescaled_batch <- plotReducedDim(rescaled,
                     dimred = "UMAP",
                     ncomponents = 2,
                     colour_by = "batch")
p1+p2

#MNN correction==========
set.seed(123)

mnn.out <- fastMNN(sce_v2, sce_v3, d=50, k=20, subset.row=HVG,
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out

#BiocSingular::RandomParam(deferred=TRUE) : complex object (similar to SSE object), contains set of parameters. Object specifying the algorithm to use for PCA.
str(BiocSingular::RandomParam(deferred=TRUE))
#class = randomParam
#Different way to run BiocSingularParam:
#IrlbaParam - calculates first PCs without calculating all
#FastAutoParam

#MNN Value: A rotation column the rowData slot, containing the rotation matrix used for the PCA.
#reconstructed matrix - once cells in single PCA space, reconstruct gene expression profile for each cell that matches the corrected PCA coordinates (see rownames(889)).
#by default it's just HVGs, but can set correct.all=TRUE for all genes.
# # in first instance, just carry out MNN on HVGs to visualise, especially with large dataset.
#Be aware of rare populations: k param depends on zie of the clusters you want. if the population of cells is 50, but you're asking for 50 NN, may take cells outs of other clusters


#Run UMAP:
mnn.out$batch <- as.factor(mnn.out$batch)

mnn.out <- runUMAP(mnn.out,
               name="UMAP",
               dimred="corrected",
               n_dimred=20)

plotReducedDim(mnn.out,
               dimred = "UMAP",
               ncomponents = 2,
               colour_by = 'batch')

#still outliers - would go back to QC and look at samples to see why they're outliers eg mt, dead etc

#Here we'll move on:

#Correction diagnostics ====
set.seed(123)
snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected")
clusters.resc <- igraph::cluster_walktrap(snn.gr)$membership
colData(mnn.out)$cluster_walktrap_SNN <- as.factor(clusters.resc)
#Set cluster to louvain
clusters.resc <- igraph::cluster_louvain(snn.gr)$membership
colData(mnn.out)$cluster_louvain <- as.factor(clusters.resc)
tab.resc <- table(Cluster=clusters.resc, Batch=mnn.out$batch)
tab.resc
colData(mnn.out)
#compare louvain and walktrap
tab.cluster <- table(Walktrap=mnn.out$cluster_walktrap_SNN, Louvain=mnn.out$cluster_louvain)
tab.cluster
#few mismatches between louvain and walktrap clustering

#heatmap of logtransformed tab.cluster
heatmap(log1p(tab.cluster))
#inspect the clustering distribution - how many cells are in both clusters or seperately? Quickly measure which cluster of louvain matches walktrap and which cell belongs to clustering method of the other.


#clustering and UMAP are independent - don't need to plot together

mnn_batch <- plotReducedDim(mnn.out,
               dimred = "UMAP",
               ncomponents = 2,
               colour_by = 'batch')

mnn_batch + rescaled_batch


#Biologically markers instead of HVGs=====

#instead of hvgs, may want to find marker genes statistically in the individual batches and then join them together. Then run MNN (on these unique markers)

#stats3 <- pairwiseWilcox(sce_v3, direction="up")
#markers3 <- getTopMarkers(stats3[[1]], stats3[[2]], n=10)

#stats4 <- pairwiseWilcox(sce_v2, direction="up")
#markers4 <- getTopMarkers(stats4[[1]], stats4[[2]], n=10)

#marker.set <- unique(unlist(c(unlist(markers3), unlist(markers4))))
#length(marker.set) # getting the total number of genes selected in this manner.

#mnn.out2 <- fastMNN(sce_v2, sce_v3, subset.row=marker.set, BSPARAM=BiocSingular::RandomParam(deferred=TRUE))


#Biologically corrected values=====
#find markers


#block (accounts for exp design)
##blocking on the batch - within each cluster - find marker adjusts for batch correct (with 1 batch reference level), adjust mean expression before calculating differential expression stats.
###Controlling for batch, testing differences between clusters.
### do batches have similar variation? - will alter p value as consequence


m.out <- findMarkers(sce, mnn.out$cluster_louvain, block=sce$Sample, direction="up", lfc=1, row.data=rowData(sce)[,3,drop=FALSE], assay.type = "logcounts")

#direction = 'up' - only keep genes that have a positive logfold change
#pairedT test by default


m.out[['2']]
 #p value = capacity of marker to separate cluster specified from all others
#logfold summary (?) look up:

## For each gene and cluster, the summary effect size is defined as the effect size from the pairwise comparison with the largest p-value. This reflects the fact that, with this approach, a gene is only as significant as its weakest DE. Again, this value is not directly used for ranking and are only reported for the sake of the user.


#Matching to symbol:

m.out[['2']]$Symbol <- rowData(sce)[rownames(m.out[['2']]), "Symbol"]
#subset DF with row (ID) with col (gene symbol)

View(as.data.frame(m.out[['2']]))



