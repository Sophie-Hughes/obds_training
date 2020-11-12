renv::restore()


library(tidyverse)
library(Seurat)
library(patchwork)
library(DropletUtils)
library(scDblFinder)

DropletUtils::read10xCounts

pbmc1k <- read10xCounts(c(pbmc1k = 'data/pbmc_1k_v3_raw/'), 
                        col.names = TRUE)

pbmc1k #33538 x 6794880  6 million barcodes in dataset and only 1000 are cells

set.seed(123)
detected_droplet <- emptyDrops(assay(pbmc1k)) #pass the matrix from the pbmc1k dataset 

detected_droplet #DataFrame with 6794880 rows and 5 columns #discard 0s
#stat with false discovery rate FDR and UMI. P value is that you're confident it's not an empty droplet. High p value for barcode means it's likely to be an empty droplet. low p value means it's likely to be a cell 

detected_droplet$rank <- rank(-detected_droplet$Total) #give 1 to the highest 

#ggplot - rank on x 

#6794880 points so will take a lot of time to run
plotdata <- as.data.frame(detected_droplet)
dim(plotdata) #6794880 x 6

#filtering
plotdata <- dplyr::filter(plotdata, Total >= 1)
dim(plotdata) #234600 x 6

#plot data, colour by FDR < 0.1
#NA = cells that have a total UMI less than 100 
#FDR is probablity cell is empty. an FDR < 0.01 is likley to be cell. FDR > 0.01 is #likely to be empty droplet 

ggplot(plotdata, aes(x = rank, y = Total, colour = FDR < 0.01)) +
  geom_point() +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") 
#blue = true cells; red likely to be empty droplets 
#lower the FDR cut off, the more cells you will exclude as empty droplets 
#If FDR < 0.01 -> more empty droplets 


pbmc1k_filter <- pbmc1k[, which(detected_droplet$FDR < 0.01)] #keep everything with an FDR less than 0.01. 
#which() tells you the positions in the vector that are TRUE, and not NA or FALSE and passes it to the function 
#NAs are all empty droplet. 

pbmc1k_filter #33538 x 1206 

#=======PART 2: doublet exercise========  

doublets <- scDblFinder(pbmc1k_filter) #1% per thousand cells 
#16 (1.3%) doublets called 
#expecting 1% of sample to be doublets in 1000 cells - 10 cells 
#
#clustering = null 

doublets
colData(doublets) #DataFrame with 1206 rows and 12 columns - adds metadata  

#check the number of singlets/doublets 
table(doublets$scDblFinder.class)
#singlet doublet 
#1190      16  

#scDblFinder.nearestClass: asks is the nearest cell another cell or artificial doublet? 
#scDblFinder.cxds_score: whats the score of the barcode - is it a cell or doublet? 

ggplot(as.data.frame(colData(doublets)), aes(x = scDblFinder.score)) +
         geom_histogram(colour = "black", fill = "white", bins = 50)
#highest scDblFinder.score is 1 (doublets)  

singlets_only <- doublets[, which(doublets$scDblFinder.class == 'singlet')]
singlets_only #dim: 33538 x 1190 (removed 16 doublets)

#summary:
#import data
#run emptyDrops
#exclude empty droplets (e.g., FDR < 0.01)
#run scDblFinder
#exclude doublets (doing this right now)


#Cellranger doesn't detect doublets - so be sure to remove these before normalising 
#load data in seurat, normalise, find variable features, scale, run PCA, run UMAP.

#After excluding empty droplets, compare how we have filtered compared to cell ranger in pbmc_filtered data 

#convert pbmc1k_filtered to seurat object - could see how many doublets that DoubletFinder - pick one of sctransform or standard method - he's overlapped two methods 
#don't use ground truth method - you won't know what's a doublet and what's a singlet 
#either assume you're 

#run 1st, 3rd and few lines from the last few sections Message by Kevin Rue
#https://github.com/chris-mcginnis-ucsf/DoubletFinder#example-code-for-real-world-applications 
#not the preferred method though! 
