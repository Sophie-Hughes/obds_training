#Perform PCA. How many principal components do you think you should keep for follow up analysis

#library import
library(ggplot2)
library(tidyverse)
library(cowplot)
library(umap)
#-------

#Exercise:  find genes associated with principle components -

#load in data
log_counts <- read.csv("data/logcounts.csv", row.names = 1)
View(log_counts)

cell_metadata <- read.csv("data/cell_metadata.csv", row.names = 1)
View(cell_metadata)

gene_metadata <- read.csv("data/gene_metadata.csv", row.names = 1)
View(gene_metadata)

#convert DF to matrix for pca analysis

log_count_matrix <- as.matrix(log_counts)
class(log_count_matrix)

#check dimension of matrix prior to PCA - know what to expect (time for PCA)
dim(log_count_matrix)

#PCA

#need to t() the data to transform it
log_count_pca <- prcomp(t(log_count_matrix), center = TRUE, scale. = TRUE)
summary(log_count_pca)

str(log_count_pca)
View(log_count_pca$x)

#plot pca

plot(log_count_pca)


#ggplot on PC1 vs PC2
ggplot(as.data.frame(log_count_pca$x), aes(x = PC1, y = PC2)) +
    geom_point()

#plot as a scree plot
screeplot(log_count_pca, npcs = 50, type='lines')

#from scree plot you might take PC10-12 (bottom of elbow, where it becomes linear)


#plot % variance
log_count_pca$sdev

#created another column with the length of the log_count_pca
#sdev is sqr root of variance so we **2 (square it)
#visualise the percentage variance

#for pc = -> could also use paste0("PC", seq_along(log_count_pca$sdev)) or
#seq(1, length(log_count_pca$sdev))

scree_pca_sdv <- data.frame(
    var = log_count_pca$sdev **2,
    pc = seq(1, length(log_count_pca$sdev))) %>%
    mutate(
        percentage_var = var/sum(var) * 100,
        cumulative_var = cumsum(percentage_var))


View(scree_pca_sdv)


#plotting percentage variance
ggplot(scree_pca_sdv[1:100, ], aes(x = pc, y = percentage_var)) +
    geom_point()


#plotting cumulative variance
ggplot(scree_pca_sdv[1:100, ], aes(x = pc, y = cumulative_var)) +
    geom_point()

#use these plots to decide how many PCs to give to UMAP

# ------------------------#
#sidenote - colour dots by experimental data


#join experimental info into pca_df DF using rownames_to_column in cell_metadata and left_join to join pca_df to cell_md

pca_df <- as.data.frame(log_count_pca$x)
View(pca_df)
pca_df <- rownames_to_column(pca_df)
View(pca_df)
cell_md <- rownames_to_column(cell_metadata)

pca_df <- pca_df %>%

    left_join(cell_md)

View(pca_df)

#plot based on experimental data to re-colour data points to see interaction with experimental conditions
ggplot(pca_df, aes(x = PC1, y = PC2, colour = Infection)) +
    geom_point()

#investigate which groups of cells are separating PCs - if we plot PC1 code for each cell for mock, STM-D, STM-L2 - can see what PC1 is separating - here it's time points for stimulated cells -
#using geom_density
#alpha = transparency

ggplot(pca_df) +
    geom_density(aes(PC1, fill = Status), color = "black", alpha = 0.5) +
    facet_grid(Time~Infection) +
    theme_cowplot()


''' could use loop to do this:

for (metadata in factor(colnames(cell_md)[-1])){
    p <- pca_df %>%
    ggplot(aes_string(x = "PC1", y = "PC2", color = metadata)) +
    geom_point()
    print(p)
    }​​
}​​ '''

#----------------------------------#
#Look in rotation to find top genes associated with PCs
View(log_count_pca$rotation)

#extract top 10 genes from PC1
topgenes_pc1 <- log_count_pca$rotation[ ,'PC1']
View(topgenes_pc1)
topgenes_pc1 #it's also added the gene names

sorted_pc1_genes <- sort(topgenes_pc1, decreasing = TRUE)
View(sorted_pc1_genes)
#ISG2 = top gene correlated with PC1 and time.

#plot expression of gene overlayed onto PCA plot (gradient left-low, right-high)

#---------------------#

#Clustering

#clustering with kmeans

kmeans <- kmeans(t(log_count_matrix), centers = 4, iter.max = 10)
kmeans
str(kmeans)

#grab clusters
#head(kmeans$cluster, 10)

#View(pca_df)
#View(kmeans$cluster)

#ordering kmeans$cluster according to pca_df$rowname
#can add 'cluster'in pca_df$cluster for the first time to add the new column
pca_df$cluster <- as.factor(kmeans$cluster[pca_df$rowname])
pca_df$cluster


ggplot(pca_df, aes(x = PC1, y = PC2, colour = cluster)) +
    geom_point()

#sum of squares
kmeans$withinss
#how distance is clustering
kmeans$betweenss

#to plot withiness, betweenss we need multiple k's using apply
candidate_k = 2:20 #DRY = don't repeat yourself
km <- sapply(candidate_k, function(i){
    km <- kmeans(t(log_count_matrix), centers = i)
    sum(kmeans$withinss)
})

km

kmeans_sum <- data_frame(sumwithinss = km, k = candidate_k)
kmeans_sum

ggplot(kmeans_sum, aes(x= k, y = sumwithinss)) +
    geom_point()

#compare cluster labels

cluster_label <- data.frame(cluster = pca_df$cluster, time = pca_df$Time)
table(cluster_label)

gg_cluster <- ggplot(pca_df, aes(x= PC1, y = PC2, colour = cluster)) +
    geom_point()

gg_time <- ggplot(pca_df, aes(x= PC1, y = PC2, colour = Time)) +
    geom_point()

gg_infection <- ggplot(pca_df, aes(x= PC1, y = PC2, colour = Infection)) +
    geom_point()

plot_grid(gg_cluster, gg_time, gg_infection)


#---------------#
#clustering using UMAP
#need to give UMAP with rotated data (x)
#-rotation is matrix of variable loading
umap_dc <- umap(log_count_pca$x)
umap_dc
#main component of object in layout - holds matrix with coordinates
umap_dc$layout
umap_coord <- as.data.frame(umap_dc$layout)
head(umap_coord)
umap_coord <- cbind(umap_coord, cell_metadata)

p1 <- ggplot(umap_coord, aes(x = V1, y = V2, colour = Time))+
    geom_point()+
    theme_cowplot()

p2<-ggplot(umap_coord, aes(x = V1, y = V2, colour = Status))+
    geom_point()+
    theme_cowplot()

p3<-ggplot(umap_coord, aes(x = V1, y = V2, colour = Infection))+
    geom_point()+
    theme_cowplot()

plot_grid(p1, p2, p3)
