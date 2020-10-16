#tidying up RNAseq data on mouse T cells.
#knockin, KO, +/-virus, 3 biological repeats (n=12)
#2 files.

library(tidyverse)
library(biomaRt)
library(pheatmap)
library(ggplot2)

#tidy count file - 3 cols
sample_table <- read_tsv("data/obds_sampletable.tsv")
count_table <- read_tsv("data/obds_countstable.tsv.gz")
View(sample_table)
View(count_table)

processed_counttable <- count_table %>%
  pivot_longer(-Geneid, names_to = "sample", values_to = "count")
View(processed_counttable)

#join with gene info to get mgi_symbol
listMarts()
ensembl <- useMart("ensembl") #connecting to specific database
datasets <- listDatasets(ensembl) #list dataset
head(datasets)
ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
View(ensembl)

filters <- listFilters(ensembl)
filters[1:5,]
View(filters)

attributes = listAttributes(ensembl)
attributes[1:5,]
View(attributes)
#matching
matching <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
      values = unique(processed_counttable$Geneid),
      mart = ensembl)
#remove duplicates
View(matching)
      
#joining tables to match geneid note different naming btw 2 tables.
processed_counttable <- processed_counttable %>% 
  left_join(matching, by = c("Geneid" = "ensembl_gene_id"))

#tidy metadata file: 1 var per col & remove species and library_layout cols
#separating by cell type, KO/KI, replicate
View(sample_table)

processed_sample_table <- sample_table %>% 
  separate(sample_title, c("gene", "genotype", "cell type", "replicates"), sep = "_") %>%
  unite("Genotype", gene, genotype, sep = "_")  %>% 
  dplyr::select(-library_layout, -read_count)
View(processed_sample_table)
#remove cols

#joining 2 tables
processed_joined <- processed_counttable %>% 
  left_join(processed_sample_table, by = c("sample" = "Sample_accession"))
View(processed_joined)

#calculate CPM (group_by() and mutate()) #log2
calculated <- processed_joined %>% 
  group_by(sample) %>% 
  mutate(total_count_per_sample = sum(count)) %>% 
  mutate(total_count_in_million = total_count_per_sample / 1000000)  %>% 
  mutate(CPM= count / total_count_in_million) %>% 
  mutate(log_CPM = log2(CPM+1))
View(calculated)

#add metadata to table w/ counts & gene info
ggplot(calculated, aes(x = sample, y = total_count_per_sample))+
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

calculated_gene_with_no_counts <- calculated %>%
  group_by(mgi_symbol) %>%
  summarise(total_count_for_each_gene = sum(count)) %>%
  filter(total_count_for_each_gene == 0) %>%
  pull(mgi_symbol)

length(unique(calculated_gene_with_no_counts))

#how many genes have a distribution of log CPM, colour by sample 
ggplot(calculated, aes(x = log_CPM, colour = sample)) +
  geom_density()

#lots of genes are very lowely expressed. So filter out genes with CPM less than, if a gene is expressed at a certain level less than 
#we've got 4 exp groups, each with 3 samples - don't want to remove genes with high expression in one group but not another. we want expression of gene in at least - we want to keep genes in one group of cells (certain time point, certain conditions) - in

calculated_filtered <- calculated %>%
  group_by(Geneid) %>%
  mutate(high_expression = sum(CPM > 0.5)) %>%
  filter(high_expression >= 3) #filter (i.e. keep) genes that have detectable expression (> 0.5 CPM) in at least 3 samples


#how many samples have CPM > 0.5 - want there to be at least 3 Trues (3 replicates)
#CPM < 0.5 is boolen vector. We want samples If all samples were less than 0.5 - we get a value of 12 as value of 12 samples in boolean vector with True (true =1) -> 12 x1 = 12 

View(calculated_filtered)

calculated$CPM < 0.5 

matching_filter <- matching %>% X
  filter(mgi_symbol == "Pakap")

#calculate the proportion of genes that were lowely expressed 

unique_gene_list <- length(unique(calculated$Geneid))  
unique_gene_list #gene before filtering

post_filtering <- length(unique(calculated_filtered$Geneid))
post_filtering #genes post filtering

calculated_proportion <- (unique_gene_list - post_filtering)/ unique_gene_list * 100  
calculated_proportion
#there's 43.5% of genes that have been filtered out. 

#make a density plot of log2(CPM+1) with filtered data 
ggplot(calculated_filtered, aes(x = log_CPM, colour = sample)) +
  geom_density() 
#can visualize we've gotten rid of lowley expressed genes on LHS of graph 

ggplot(calculated_filtered, aes(x = log_CPM, colour = sample)) +
  geom_histogram() +
  facet_wrap(~sample) #facet by 1 variable, facet_grid facets on x AND y

#split by sample/plot 

#-----------------------#

#Plot CD4 and CD8 expression for all samples:
#CD4/CD8 on x axis with expression on y, 1 point per sample to see diff in expression across samples 

filtered_t_genes <- calculated_filtered %>%
  filter(mgi_symbol == "Cd4" | mgi_symbol == "Cd8a") 
  #filter(mgi_symbol %in% c("Cd4, Cd8a")) - alternative method to filter CD4/CD8a

View(filtered_t_genes)

ggplot(filtered_t_genes, aes(x = mgi_symbol, y= log_CPM, fill = replicates)) +
  geom_col(position = "dodge") +
  facet_grid(Genotype ~ `cell type`)  # use `` when you have spaces 

#filter Erg2/3 and plot 

filtered_egr_genes <- calculated_filtered %>% 
  filter(mgi_symbol == "Egr2" | mgi_symbol == "Egr3") 

View(filtered_egr_genes)
View(calculated_filtered)

ggplot(filtered_egr_genes, aes(x = mgi_symbol, y= log_CPM, fill = replicates)) +
  geom_col(position = "dodge") +
  facet_grid(Genotype ~ `cell type`) 

#Choose 8 biologically relevant genes and plot heatmap 
relevant_genes <- calculated_filtered %>%
  group_by(Geneid) %>%
  summarise(variance = var(log_CPM)) %>% #summarise end up 1 value/gene
  arrange(desc(variance)) %>%   #tidyverse sort function 
  slice_head(n = 8)

View(relevant_genes)

#keep 8 genes: 
relevant_genes_count <- calculated_filtered %>%
  filter(Geneid %in% relevant_genes$Geneid) 

unique(relevant_genes_count$Geneid)
         
#heatmap 
#convert to table to swap direction of col and row (vs count_table)
#filter() for row
#select() for col 

count_table_heatmap <- relevant_genes_count %>%
  ungroup() %>% #need to ungroup, tidyverse remembers groups until you tell it to forget
  dplyr::select(mgi_symbol, sample, count) %>%
  pivot_wider(names_from = sample, values_from = count) %>%
  column_to_rownames(var = "mgi_symbol") #take id col and make row as need a matrix for heatmap (all numeric)

View(count_table_heatmap)
#pheatmap should need a matrix, but is clever and converts DF to matrix. Otherwise use pheatmap(as.matrix(count_table_heatmap)) 

pheatmap(count_table_heatmap)
#first 5 genes drawfs other gene expression 

pheatmap(count_table_heatmap, scale = 'row') 
#scale = 'row' allows us to see variance better across samples by scale-ing across genes 
#scale = remove mean and divide by standard deviation 
#gene is expressed 1.5 standard deviations 
#genes are independent of each other 
#z-score relevant to mean of gene 
#z-score is number of sd away of mean 
#zscore = 0 -> same value as mean 
#when its scaled we can compare samples, but can't compare multiple genes within a sample - read rows but can't read columns 

#here we can say IL-21 has high expression in first 3 samples. but better to do scalling by gene not by sample 
