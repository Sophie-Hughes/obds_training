#restore environment 
renv::restore()

#library import 
library(EnsDb.Mmusculus.v79) 
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
#read in data
data <- read.csv("data/resApeg_padj.csv")
dim(data)


#we want to generate 2 seperate gene table: 
#upregualted: with padj < 0.05 and log2FC > 1 
#downregulated:with padj > 0.05 and log2FC  1  

#Get GENEID and ENTEREZ ID 

rownames(data)<-data$GENEID
edb <- EnsDb.Mmusculus.v79
columns(edb)
gene_id <- select(edb, keys = rownames(data), columns = c("GENEID", "ENTREZID"), keytype = "GENEID")
View(gene_id)

data <- data %>%
    left_join(gene_id)

#----------------#

upregulated_genes <- data[data$padj < 0.05 & data$log2FoldChange > 1,]  #Keep all cols
downregulated_genes <- data[data$padj < 0.05 & data$log2FoldChange < -1,] 

#filter out NA ENTEREZID values:
upregulated_genes_withoutNA <- upregulated_genes %>% 
    filter(!is.na(ENTREZID))

downregulated_genes_withoutNA <- downregulated_genes %>% 
    filter(!is.na(ENTREZID))


## feature 1: numeric vector
genelist_up <- upregulated_genes_withoutNA[,3]

## feature 2: named vector
names(genelist_up) <- Upregulated_genes_withouNA$ENTREZID

## feature 3: decreasing order
genelist_up <- sort(genelist_up, decreasing = TRUE)

#----------------#

data_filtered <- data %>% 
    filter(!is.na(ENTREZID))

data_filtered_list <- data_filtered[,3]
names(data_filtered_list) <- data_filtered$ENTREZID
data_filtered_list <- sort(data_filtered_list, decreasing = TRUE)


#choose method - GO over-representation test (5.3)

#feed in vector of entrez geneid into genes = in enrichGO function. We've defined it outside the enrichGO argument 
names(genelist_up) <-upregulated_genes_withoutNA$ENTREZID
table(names(genelist_up)) 
head(sort(table(names(genelist_up), useNA ='always'), decreasing = TRUE)) #1 gene is duplicated 


GO <- enrichGO(gene = names(genelist_up),
               universe = names(data_filtered_list),
               OrgDb = 'org.Mm.eg.db'
               ) 


#plotting 
barplot(height = GO, x = 'GeneRatio', color = "p.adjust", showCategory = 10)

#gene_concept_mapping
matching_go <- setReadable(GO, "org.Mm.eg.db",'ENTREZID')
cnetplot(x = matching_go, layout = 'nicely')

#Enrichment map 
emapplot(GO, showCategory = 20, pie_scale = 1.5, layout = "nicely")

#dotplot 
plot1 <- dotplot(GO, showCategory = 20) 

plot1 + ggtitle('dotplot for GO') +
    theme(axis.text.x = element_text(size = 7))

#give plots back as ggplots so we can use ggplot functions 

#KEG analysis: 

KEGG <- enrichKEGG(gene = names(genelist_up),
                 organism = 'mmu',
                 pvalueCutoff = 0.05,
                 universe = names(data_filtered_list))


# use https://www.genome.jp/kegg/catalog/org_list.html to find organism ID
#use default padj = BH (so don't need to specify)

class(KEGG) #enrichResult, dose

KEGG_result <- KEGG@result
View(KEGG_result)
sum(duplicated(names(genelist_up)))

#plotting 
barplot(height = KEGG, x = 'GeneRatio', color = "p.adjust", showCategory = 10)

#gene_concept mapping
matching <- setReadable(KEGG, "org.Mm.eg.db",'ENTREZID')
cnetplot(x = matching, layout = 'nicely')

#Enrichment map 
emapplot(KEGG, showCategory = 20, pie_scale = 1.5, layout = "nicely")

#-------------------#

#had error message when running gseGO to remove the duplicates: 
#remove duplicated names 
gse_clear <- data_filtered_list[!duplicated(names(data_filtered_list))] 

length(gsa_clear)
length(data_filtered_list)

#GO gene set enrichment analysis
#warning message with tied genes is fine. 

go_enrich <- gseGO(geneList = gse_clear,
                   OrgDb = org.Mm.eg.db,
                   ont = 'ALL',
                   minGSSize = 100,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   verbose = TRUE)

gsa$Description[1:3] 
                
#plot the enrichment 

ridgeplot(go_enrich)

enrichplot::gseaplot2(go_enrich, geneSetID = c("GO:0045087","GO:0002768"))

terms <- go_enrich$Description[1:3] enrichplot::pmcplot(terms, 2010:2020)
