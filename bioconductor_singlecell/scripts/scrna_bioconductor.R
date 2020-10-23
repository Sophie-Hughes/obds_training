#restore environment
renv::restore()

#import library
library(tidyverse)
library(scran)
library(scater)
library(DropletUtils)
library(DelayedMatrixStats)
library(DelayedArray)
library(Matrix)
library(sparseMatrixStats)
library(iSEE)
library(biomaRt)


# Give folder with the barcodes.tsv.gz, matrix.mtx.gz, features.tsv.gz
sc_data <- read10xCounts(c(filtered = "data/filtered_feature_bc_matrix"), col.names = TRUE)


sc_data #SingleCellExperiment class
str(sc_data)
slotNames(sc_data) #tells the number of slots available

#look at the object:
colData(sc_data)
colnames(colData(sc_data)) #Sample and Barcode
rowData(sc_data)
colnames(rowData(sc_data)) #ID, Symbol, Type
assayNames(sc_data) #Counts
metadata(sc_data) #access metadata
assay(sc_data) #provides sparse matrix
# # when multiple assays are laoded then to select which assay to use: assay(sc_data, 'counts') or assay(sc_data)$counts


#sample is path to filtered matrix. Useful if you have multiple samples in a single object. Can read in sample names

#Compute and visualise quality control metrics; use scater. Would you remove any cell?

#Barcode ranks needs a sparse matrix
barcode <- barcodeRanks(assay(sc_data)) #note we're working on a pre-filtered dataset

ggplot(data.frame(barcode), aes(x=rank, y=total, col=fitted)) +
    geom_point()+
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(trans = 'log10')

#using scater: per cell QC:
per_cell <- perCellQCMetrics(sc_data)

#argument: percent_top (% UMI accounted for by the top n genes, by default it looks at range) Looking for cell with diverse transcriptome- 100% of reads come from as many genes as possible. (Specify name of count matrix, 'counts' by default - many need to change this argument depending on what's in object eg '10X counts')

#top 50 genes histogram
ggplot(as.data.frame(sc_data), aes(x = percent_top_50)) +
    geom_histogram(bins = 50, colour = 'black')

#for around 30 cells, 100% of of the cells account for the top 50 genes
#top 50 genes account for 30-40% of library size.
#any cell were top 50 genes account for 60% of library = poor quality because of low counts and small transcriptome (?) but may want to interogate this with further scatter:
ggplot(data.frame(sc_data), aes(x = percent_top_50, y = total)) +
    geom_point() +
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(trans = 'log10')

#above 60% of library has: low total count and high % of top 50 genes - but these genes are most highly expressed in the library.
#keep visualising and make a decision about where to cut your threshold. We'll cut off at 60%
#could also colour by high mt_metric to see whether these high reads may be dying cells with high mt (mitochondria)




ggplot(data.frame(sc_data))+
    geom_col(aes(y=sum,x=rownames(data.frame(sc_data))), bins=100, colour="black") +
    geom_col(aes(y=detected, x=rownames(data.frame(sc_data))), bins=100, colour="red")


#plot all percent top genes in one go:

sc_data %>% as.data.frame %>%
    pivot_longer(-c(sum, detected, total)) %>%
    mutate(name = factor(name, levels = unique(name))) %>%
    ggplot(aes(x=value, y = total)) +
    geom_point() +
    facet_wrap(~name)

#Another per cell QC
sc_data <- addPerCellQC(sc_data)
sc_data #single cell experiment object
colData(sc_data) #added all col from DF into field.

ggplot(as.data.frame(colData(sc_data)), aes(x = percent_top_50)) +
    geom_histogram(colour = 'black', fill = 'white', bins = 50)


#perfeature QC
feature_qc <- perFeatureQCMetrics(sc_data)
feature_qc #mean count/feature & detected (fraction of samples that genes have been detected)
sc_data <- addPerFeatureQC(sc_data) #'stuffing' into same object
rowData(sc_data) #add mean & detected (n=5cols)

#we are NOT filtering yet!
#neutrophils are difficult to sequence, but this may be because they're filtered out
# # so don't filter too soon, be sure!

#LogNorm===============

#Convert the counts into normalized expression values to eliminate cell-specific biases (e.g. in capture efficiency) use scater and scran

#preprocessing: log normalised expression values
sc_data <- logNormCounts(sc_data)
sc_data


#lognormcount calculating log normalised for each cell then using log calculated matrix on whole library - adjusting for library size and then log transforming. Equivalent of counts per million in bulk RNAseq Uusally counts/10,000 in scRNA

range(assay(sc_data, 'logcounts')) #0.00000 - 12.68275

rowMeans(assay(sc_data, 'logcounts')) #need to specify 'assay' otherwise will pick the 1st 'count'
sparseMatrixStats::rowVars(assay(sc_data, 'logcounts')) #

mean_var <- data.frame(mean = rowMeans(assay(sc_data, 'counts') + 1),
           variance = sparseMatrixStats::rowVars(assay(sc_data, 'counts') + 1)
           )

#mean as 1st col:
#mean on x, var on y
ggplot(data.frame(mean_var), aes(x= mean, y= variance))+
    geom_point() #post-normalisation

ggplot(data.frame(mean_var), aes(x= mean, y = variance))+
    geom_point() #lower mean and var pre-normalisation

#Filtering================
#Select features for downstream analyses, e.g. highly variable genes; use scran.

decomposed_log <- modelGeneVar(sc_data) #modelling var of log-exp for each gene
decomposed_log #technical/biological component of the var
#apply DR to compact data nad further reduce noise and scatter

plot(decomposed_log$mean, decomposed_log$total)
curve(metadata(decomposed_log)$trend(x), add=TRUE, col="dodgerblue")
#each dot = gene; line = trend of total variance -> trend is where most of genes are
# # Assume technical variance is where majority of genes are
# # genes taht show more variance than typical (top 3 at logtotal 8) they are ariable in data set as a whole - want to find these in unsupervised way. these genes are driving heterogeneity but we want to grab enough genes so:
# # we want to grab genes above blue line (biological variance); biological varaince (above blue line) is beyond technical varaince (blue line).

#Feature Selection:
#define HVGs (above blue line)
hvgs <- getTopHVGs(decomposed_log,
           var.field = 'bio',
           n = 1000)

#taking biological variance which is it's own column 'bio' which we can select

rowData(sc_data)$HVG <- rownames(sc_data) %in% hvgs


#dimensionality========
#Apply dimensionality reduction to compact the data and further reduce noise; scater

#set seed for PCA reproducibility
set.seed(123)

sc_data <- runPCA(sc_data, name = 'PCA', ntop = Inf, subset_row = hvgs)
#use all HVGs with Inf

#runPCA returns PCA result as rowname in singlecellexperiment object
#calculatePCA returns the PCA directly and is run first. PCA is only returning first 50 PCs.

reducedDims(sc_data)

str(reducedDims(sc_data)$PCA) #give PCA component from object
#matrix with each col per Principle component, each row is a cell ?

plotReducedDim(sc_data, dimred = 'PCA', colour_by = 'ENSG00000090382')
#top 2 PCs can explain more of the variance

#Hacky version to plot scre elbow plot
percent.var <- attr(reducedDim(sc_data), "percentVar")
plot(percent.var, xlab="PC", ylab="Variance explained (%)")



#runUMAP========
sc_data <- runUMAP(sc_data, name = 'UMAP', dimred = 'PCA', n_dimred = 30)

#specify the dimred you've already computed - 'PCA'
#Can specify which dimensions you want it to us with n_dimred

plotReducedDim(sc_data, dimred = "UMAP", colour_by = 'ENSG00000090382')
colData(sc_data)

#print(head(hvgs)) to find gene names
#can colour by gene name on PCA and UMAP too

#cluster==========
#Cluster cells; use scran


#Quick cluster - is less precise:
colData(sc_data)$quickCluster <- quickCluster(sc_data, assay.type="logcounts")

plotReducedDim(sc_data,
               dimred = "UMAP",
               ncomponents = 2,
               colour_by="quickCluster")

SNNGraph <- buildSNNGraph(sc_data, use.dimred = 'PCA')
colData(sc_data)[["cluster_louvain"]] <- factor(igraph::cluster_louvain(SNNGraph)$membership)

colData(sc_data)[["cluster_louvain"]]

plotReducedDim(sc_data,
               dimred = 'UMAP',
               colour_by = 'cluster_louvain')

#identify_markers========
#identify markers per cluster using scran

markers <-findMarkers(sc_data, groups = sc_data$cluster_louvain)
markers[[1]] #take 1st element of list to compare this cluster with other clusters

#comparing this gene to other clusters

iSEE(sc_data)

#Matching to symbol:
rowData(sc_data)
rownames(sc_data) <- uniquifyFeatureNames(rownames(sc_data), rowData(sc_data)$Symbol) #provide ensemble gene ID and gene symbol

#can now see geneID
iSEE(sc_data)

