
#library imports

library(tidyverse)
library(biomaRt)
library(pheatmap)
library(DESeq2)
library(cowplot)
library(patchwork)
library(EnsDb.Mmusculus.v79) 

#renv::install('cowplot')
#renv::install('patchwork')


#read in data and metadata 

sample_table <- read_tsv("data/obds_sampletable.tsv")
counts_table <- read_tsv("data/obds_countstable.tsv.gz")
View(sample_table)
View(counts_table)

#Convert the counts table (obds_countstable.tsv.gz) and the sample information table (obds_sampletable.tsv) into a suitable format for generating a DESeqDataSet object
class(counts_table)
counts_table <- column_to_rownames(counts_table, "Geneid")
counts_table <- as.matrix(counts_table)

class(sample_table)
sample_table <- column_to_rownames(sample_table, "Sample_accession")
class(sample_table)


#do the colnames of count_table == to rownames of sample table 
table(colnames(counts_table)==rownames(sample_table))

#Set Egr2/3 DKO CD8 cells as the reference level
#Separate sample_title column

sample_table <- sample_table %>%
    separate(sample_title, c("egr_locus", "genotype", "cell_type", "replicate"), sep = "_") %>%
    unite(col = "condition", egr_locus, genotype, cell_type, sep = "_") %>%
    dplyr::select(-c(species, library_layout)) %>%
    mutate(condition = factor(condition, levels = c("Egr2/3_DKO_CD8", "Egr2/3_DKO_CD4", "Egr2_Kin_CD4", "Egr2_Kin_CD8"))) 
sample_table
levels(sample_table$condition)

#Generate a DESeqDataSet object named dds
#Make DDS from matrix 
dds <- DESeqDataSetFromMatrix(counts_table,
                       sample_table, 
                       ~ condition)

colData(dds)
rowRanges(dds)


#Access the design formula, counts matrix and sample information from dds
design(dds) #formula 
counts(dds) #count matrix - equivalent to assays(dds)$counts


#Calculate the size factors for each sample
dds <- estimateSizeFactors(dds)
sizeFactors(dds) #counts divided by that value 

#sizeFactorsdf <- data.frame(sizeFactors(dds))

sizeFactorsdf <- data.frame(sample = names(sizeFactors(dds)),
                           size_factor = sizeFactors(dds),
                           sample_group = colData(dds)$condition)

#Generate a bar plot of the size factors for each sample, coloured by condition/group
sizeFactorsdf %>%
    ggplot(aes(y=size_factor, x= sample, fill =sample_group)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
#variability in library size that isn't affected by condition 

#Obtain dispersion estimates for each gene

dds <- estimateDispersions(dds) #only change params they recommend 
dispersions(dds)
plotDispEsts(dds) #can see dispersion affected by mean counts 
#fit dispersion towards line apart from genes above trend-line - they keep gene calculated dispersion 

#Perform the Wald test
dds <- nbinomWaldTest(dds)

#Use the DESeq() function to perform steps (estimate size factor, estimate dispersion, wald test )
#remember to run dds function again 

dds <- DESeq(dds) 

#Access the coefficients of the NB GLM
head(coef(dds), 3) 
# - NA may represent independently filtered out genes, or several other reasons #outlier - may be removed if only a few replicates 
#intercept - double KO Egr2/3 
#set baseline factor condition as Egr2/3 double kO CD8
#for each gene it's the log expression of baseline 
#change of expression in CD4


#Access the results table for the comparison between CD8+ and CD4+ T cells from Egr2/3 DKO mice
res <- results(dds, contrast = c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'))
res <- results(dds, name = 'condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8') 

           
#results going into a different type of object 
res_df <- as.data.frame(res)
View(res_df) 
#basemean - average condition across all samples 
#padjust by BH

#Plot a histogram of the raw and BH-adjusted p-values – do they look as expected? 

padj_plot <- ggplot(res_df, aes(x = padj)) +
    geom_histogram() 

pvalue <- ggplot(res_df, aes(x = pvalue)) +
    geom_histogram()


plot_grid(pvalue, padj_plot) #plot together with cowplot 

padj_plot / pvalue #plot together with patchwork 


#Generate an MA plot of the log2 FC values for all genes

plotMA(res, ) #comparing CD4 to CD8 DKO 

#Shrink the log2 FC values using the normal, apeglm and ashr methods 

#contrast and coef both select the same thing 
lfcShrink(res, coef = 'condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8')
resNorm <- lfcShrink(dds, contrast = c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'), type = 'normal')

resAsh<-  lfcShrink(dds, contrast = c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'), type = 'ashr')

resApeg<-  lfcShrink(dds, coef = 'condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8', type = 'apeglm')

#compare methods with MA plot
plotMA(resAsh)
resAshplot <- recordPlot()

plotMA(resNorm) 
resNormplot <- recordPlot() 

plotMA(resApeg) #typically use this as it's default 
resApegplot <- recordPlot()

#plot_grid(resAshplot, resNormplot, resApegplot, nrow = 3) plot_grid won't work here 
#as plotMA is baseR, can store using recordPlot

# Generate a results table (one shrinkage method) containing mgi symbols 

resApeg_table <- as.data.frame(resApeg)
View(resApeg_table)

#Use the EnsDb.Mmusculus.v79 package 
edb <- EnsDb.Mmusculus.v79
edb
columns(edb)
gene_id <- select(edb, keys = rownames(resApeg_table), columns = c("GENEID", "SYMBOL"), keytype = 'GENEID') 
View(gene_id)

#remove duplicated gene IDs 
gene_id$GENEID[duplicated(gene_id$GENEID)] #access column and subset by boolean vector
gene_id$SYMBOL[duplicated(gene_id$SYMBOL)]

#Remove all genes with a padj of NA 

duplicated_genes <- dplyr::filter(gene_id, SYMBOL %in% gene_id$SYMBOL[duplicated(gene_id$SYMBOL)])
View(duplicated_genes)

#move ensemble ID from rownames to a separate column 
resApeg_table <- rownames_to_column(resApeg_table, "GENEID") %>%
    left_join(gene_id)
View(gene_id)

#remove duplicate symbols 
resApeg_table[is.na(resApeg_table$SYMBOL),]
sum(is.na(resApeg_table$SYMBOL))


#Remove all genes with a padj of NA
resApeg_padj <- dplyr::filter(resApeg_table, !is.na(padj))
sum(is.na(resApeg_padj$padj))

#resApeg_table <- head(resApeg_table[!is.na(resApeg_table$padj),]) #alternative method

#Write the results table to a CSV file
write.csv(resApeg_padj, file = "results/resApeg_padj.csv", quote = FALSE, row.names = FALSE)

#Filter the results table for padj < 0.05, and write to a CSV file

filtered_table <- dplyr::filter(resApeg_padj, padj < 0.05) %>%
    dplyr::filter(abs(log2FoldChange) > 1) 
dim(filtered_table)

write.csv(filtered_table, file = "results/filtered_table.csv", quote = FALSE, row.names = FALSE)

#Generate VST and rlog transformed counts: 
#see relationship between mean and var - we want to find genes not because they're highly variable but because there's a diff between conditions. So we use VST to get rid of the relationship between mean and var 

vsd <- vst(dds, blind = FALSE)
#generally use blind = false - knows which condition your sample is in. 

rld <- rlog(dds, blind = FALSE) 

#Plot the relationship between the mean expression (X-axis) and the sd of all genes (y-axis) 

vsn::meanSdPlot(assay(vsd))
#trend line is flat, var is different between low and high 

vsn::meanSdPlot(assay(rld))
#rlog looks odd, so we'll take vst counts 

#Generate a PCA plot either using all genes, or top 500 most variable genes 
#we can use deseq2 PCA here 

#plot all genes 
plotPCA(vsd, ntop = nrow(vsd))

#clear separation 

#plot 500 genes
plotPCA(vsd)
#we want to colour by condition which is default so don't need to specify intgroup = here 
#ntop = 500 is the default 
#top500 has even more var - max diff between groups and min differnce between groups 
#tend to plot with all genes 

#Generate a heatmap of the top 20 (by shrunken FC) differentially-expressed genes – label samples by condition and genes by mgi symbol 

#to get top 20 we need to sort by Log2FoldChange
top_20 <- filtered_table[order(-abs(filtered_table$log2FoldChange)),] %>%
    slice_head(n=20) %>%
    pull(GENEID, SYMBOL)
#order rows but col stay same, order descending using -abs 

vst_top20 <- as.data.frame(assay(vsd)) %>% 
    dplyr::filter(rownames(.) %in% top_20)

pheatmap(vst_top20, scale = "row")

#Volcano plot

df_for_volcanoplot <- resApeg_padj %>%
    mutate(y=-log10(padj)) %>%
    mutate(highlight = (padj < 0.5 & abs(log2FoldChange)>1))

colour <- c('black', 'red')
highlight <- with(df_for_volcanoplot, factor(padj < 0.5 & abs(log2FoldChange)>1))
colour <- colour[highlight]

ggplot(df_for_volcanoplot, aes(x=log2FoldChange, y=y, colour=highlight)) + 
    geom_point(colour=colour)

top_genes <- names(top_20)[!is.na(names(top_20))]
df_for_volcanoplot_subset <- df_for_volcanoplot %>%
    dplyr::filter(SYMBOL %in% top_genes)

df_for_volcanoplot <- df_for_volcanoplot %>% 
    mutate(y=-log10(padj)) %>%
    mutate(highlight = (padj < 0.5 & abs(log2FoldChange)>1))

ggplot(df_for_volcanoplot, aes(x=log2FoldChange, y=y)) +
    geom_point(colour=colour) +
    geom_label(aes(label = SYMBOL), data = df_for_volcanoplot_subset)
