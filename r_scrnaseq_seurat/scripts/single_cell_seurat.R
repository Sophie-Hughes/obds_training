renv::restore()

library(tidyverse)
library(Seurat)
library(patchwork)
library(RcppArmadillo)
library(Seurat)



# features <- read_tsv("/Users/andresnoe/obds_sep20/working_directory/filtered_feature_bc_matrix/features.tsv.gz")
# matrix <- read_tsv("/Users/andresnoe/obds_sep20/working_directory/filtered_feature_bc_matrix/matrix.mtx.gz", comment = "%")
# barcodes <- read_tsv("/Users/andresnoe/obds_sep20/working_directory/filtered_feature_bc_matrix/barcodes.tsv.gz")
tail(features)
head(matrix)
head(barcodes)

# Give folder with the barcodes.tsv.gz, matrix.mtx.gz, features.tsv.gz
pbmc.data <- Read10X("data/")

pbmc.data
View(pbmc.data)
class(pbmc.data) 
summary(pbmc.data)


#Create Seurat object 
#create assay object with ADT
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`,
                           project = "pbmc3k")
pbmc[["ADT"]] <- CreateAssayObject(counts = pbmc.data$`Antibody Capture`)

Assays(pbmc)
DefaultAssay(pbmc) #RNA is the default already, don't need to change it yet

# Initial quality control

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#store experimental metadata (which contains TCR, and BCR)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") #grep for any genes that contain "MT-".
# Can also do for any other genes that contain a string in each cell
pbmc[["percent.mt"]]

pbmc[["percent.rib"]] <- PercentageFeatureSet(pbmc, pattern = "^RPS|^RPL")
pbmc[["percent.rib"]]

View(pbmc[[]])

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA")) +
    geom_hline(yintercept=1250)

VlnPlot(pbmc, features = c("nCount_RNA"))

VlnPlot(pbmc, features = c("percent.mt")) +
    geom_hline(yintercept=14)

VlnPlot(pbmc, features = c("percent.rib"))

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC plot of nCount_RNA by nFeature_RNA coloured by percent.mt (use ggplot2)
ggplot(pbmc[[]], aes(x=nCount_RNA, y=nFeature_RNA, col=percent.mt))+
    geom_point()
# Suggestion: if just looking at this plot, can cut off at nFeature_RNA = 6000
# But this is an iterative process

# Cutoffs:
# # 1000 for nFeature_RNA (total number of RNAS)
# # 12.5 for percent.mt
# # Leave nCount_RNA for now (total number of counts) 
# # Leave percent.rib for now (but some people say that if it's above 50% in PBMCs then they should be filtered out, but be careful)

# WILL NOT FILTER YET - get up to clustering step and look for cluster of dead cells

#Remeber that SCTRansform performs these three functions:
#normalise data 
#find variablesfeatures
#scale data 


#pbmc <- NormalizeData(pbmc, normalization.method = 'LogNormalize', scale.factor = 100000)

pbmc <- SCTransform(pbmc, 
                    assay = 'RNA', 
                    seed.use = 1448145,
                    verbose = TRUE)
#this produces a new assya which has all the sc transfrom stuff
#Look at the documentation of SCTransform for defaults, particularly, 
#new.assay.name = 'SCT'
#variable.feature.n = 3000
#vars.to.regress = NULL (can regress out percent.mt) 
#return.only.var.genes = TRUE (only returns genes that are highly varibale genes)
# # might be important if gene OUTSIDE of this highly varible gene list is different experiment 
#don't regress out batch affects here (use integration to deal with them), here we're regressing out variables that have an effect but not a huge effect. 


#Warning: Don't worry about this one: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached -> seems to be fine

VariableFeature(pbmc) #find the highly variable genes and extract them as characters
#identifies which genes have been pulled by SCTransform 


#perform linear dimensional reduction 
DefaultAssay(pbmc) #SCtransform automatically changes default assay for SCT 

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

DimPlot(pbmc, reduction = "pca", group.by = 'percent.mt') +
    NoLegend()

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#print the PC loadings: genes that most contribute to pC 

VizDimLoadings(pbmc,dims = 1:3, reduction = 'pca')
#plots above print statement essentially 
#compare across genes 
#argument balanced = default is False 
#balanced = false return an equal number of genes with + and - scores. If false (default), returns the top genes ranked by the scores absolute values.

DimHeatmap(pbmc, dims = 1, cells = 2000, balanced = TRUE)
#by PC14 it's gone - the more noise there is, the less interesting that PC is likely to be
#Dimheatmap can give an indication of which PCs we should include in the dim red (similar to Scree plot)
#each column is a cell, each colour is a gene 

#Determine the 'dimensionality' of the dataset 

#JackStraw cannot be run on SCTransform-normalized data.
#pbmc <- JackStraw(pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

ElbowPlot(pbmc, ndims = 50)
#From this elbow - we'd take 25 PCs
#using SCTransform can include a few more PCs than otherwise, so err on side of increasing 
#in this case we will chose between 20-30 (25) 

#Cluster the cells 
#input to cluster is PCA space, build graph of distance in PCA space, 

pbmc <- FindNeighbors(pbmc, dims = 1:25, k.param = 20)
#dims = 25 PCs chosen from elbow plot 
#k.param play around with this, (defautlt is 20), depending on size of clusters 
#annoy.metric: distance metric fro anoy. options include euclidean, cosine, manhattan, and hammimg. 
# # ecluidean is better for immuen cells 
# # just a distance metric for measuring distance between cells 
# # annoy = approx nearest neighbors (oh yeah) is another R package 

pbmc <- FindClusters(pbmc, resolution = 0.5) 
#run different resolutions 
#We'll use Clustree (mulptiple resolution clustering). Resolution of clustering is how much you want to splot up data. higher resolution = more clusters. Run multiple resolutions to see if clusters are biologically meaningful, then choose best resolution. No best answer - just try to choose best resolution that fits biology 
#From vignette: we find that setting this parameter [resolution] between 0.4 - 1,2 typically returns good results dor single cell datasets of around 3K cells
# # affected by the nUMBER OF CELLs you have 
# # usually use clsuter tree to determine more definitely, but this is just a start 
pbmc <- FindClusters(pbmc, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.3, 1.5, 2))  
View(pbmc[[]])
Idents(pbmc) <- "SCT_snn_res.0.8"
    
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


#Run non-linear dimensional reduction (UMAP/tSNE) 

pbmc <- RunUMAP(pbmc, dims = 1:10)

#note that you can set 'label =TRUE' or use the LabelClusters function to help label 
#individual clusters
DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = nCount_RNA > 11000)) 

DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = nFeature_RNA < 1000)) 

DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = percent.mt < 12.5))

#thresholds look good as deads cells are clustering together. 
#We can filter our data to get rid of dead cells 

ggplot(pbmc[[]], aes(x=nFeature_RNA))+
    geom_histogram(bins = 100)
#from this ggplot histogram we've removed very top cells (5500) in the below subset 

pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 1000 & nFeature_RNA < 5500 & percent.mt < 12.5)

#once you have your thresholds - rerun everything from SCTransform 

#Filtering=========
#below I have copied the code from the above (in section titlted "not yet filtering")
pbmc <- SCTransform(pbmc, 
                    assay = 'RNA', 
                    seed.use = 1448145,
                    verbose = TRUE)



#=============

pbmc <- RunUMAP(pbmc, dims = 1:10) 

DimPlot(pbmc, reduction = 'umap')
#B cells bottom left?
#top monocytes
       
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))         

#post-filtering steps========= 
#lets find out how many clusters we should have 
library(clustree)
clustree(pbmc)
#but we would usually go through and look at markers of the clusters
#lets choose a resolution of 0.8

Idents(pbmc) <- "SCT_snn_res.0.8"

#subsetted PBMC = 3988 

#Add the protein expression levels to the Seurat object

#protein data also needs to be normalised and scaled like genes, but differently. we don't use SCTransform, but we're using method advised by Seurat
# # there's a new version of normalising on Suerat 

#normalise by CLR method
#scaled by -mean *SD for each antibody 
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, assay = "ADT")

#reset defualt to ADT to see what's in the object 
DefaultAssay(pbmc) <- 'RNA'
DefaultAssay(pbmc) <- 'ADT'
rownames(pbmc)

FeaturePlot(pbmc, features = c("CD3-TotalSeqB", 'CD3E',"CD4-TotalSeqB", 'CD4', "CD34-TotalSeqB", 'CD34', "CD8a-TotalSeqB", 'CD8a'), max.cutoff = "q95", ncol = 3)

RidgePlot(pbmc, features = c("CD3-TotalSeqB", "CD45RA-TotalSeqB", "CD20-TotalSeqB", "PD-1-TotalSeqB"), ncol = 2)


#Cluster directly on protein levels - PCA directly onto protein 

#perform PCA directly onto proteins
#ensure DefaultAssay set to ADT 
rownames(pbmc)
IgG <- c("IgG1-control-TotalSeqB", "IgG2a-control-TotalSeqB", "IgG2b-control-TotalSeqB")

pbmc <- RunPCA(pbmc, features = rownames(pbmc)[!rownames(pbmc) %in% IgG], npcs = 20, reduction.name = "pca_adt", reduction.key = "pca_adt_", verbose = FALSE)


DimPlot(pbmc, reduction = "pca_adt")

#run elbow plot to get dimensions 
ElbowPlot(pbmc, ndims = 15)


pbmc <- RunUMAP(pbmc, dims = 1:10, reduction = "pca_adt", reduction.key = "adtUMAP_", reduction.name = 'umap')

DimPlot(pbmc, reduction = "umap_adt") #coloured by RNA expression 


pbmc <- FindNeighbors(pbmc, features = rownames(pbmc)[!rownames(pbmc) %in% IgG], dims = NULL)
#features only used when dims = NULL, build graph on features when you have a small number of proteins. you would'nt do this for genes as there would be too many. Try with PCA dimensions (genes) or give all the proteins (features). 

pbmc <- FindClusters(pbmc, resolution = c(0.1, 0.2, 0.3, 0.4), graph.name = "ADT_snn")

#we can compare the RNA and protein clustering, and use this to annotate the protien clustering 
#we could also of course use FindMarkers

Idents(pbmc) <- 'ADT_snn.0.2'

#compare RNA and ADT clusters to look at how well they correlate 
clustering.table <- table(Idents(pbmc), pbmc@meta.data$SCT_snn_res.0.8)
clustering.table
umap_adtClusters <- Dimplot(pbmc, reduction = )
