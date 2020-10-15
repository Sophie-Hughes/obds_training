library(tidyverse)
library(gridExtra)
library(cowplot)
library(patchwork)
library(scales)
library(biomaRt)
data("mtcars")

str(mtcars)

#Draw a bar plot showing the number of cars with 3, 4 or 5 forward gears 
#Change plot margins to 6, 6, 5, 5
#Add x and y axis labels
#Add a main title spread across two lines Change the colour and width of the bars Add a #horizontal line at y = 6
par()$lwd
par(mar=c(6,6,5,5), lwd= 5, mfrow=c(1,2))


table(mtcars$gear)
barplot(table(mtcars$gear),
        xlab = "Number of gears",
        ylab = "Number of cars",
        main = "A main title\nspread across two lines",
        col = "red",
        border = "blue")
abline(h=6, lwd=7)

#\n to spread over two lines 

#Generate a scatter plot of mpg vs. hp coloured by gear values

#making vectors of colours, named by gears. 
colours <- c("red", "green", "blue")
names(colours) <- c(3,4,5)
colours
gear_column <- as.character(mtcars$gear)
colours[gear_column]
gear_column

par(mar=c(6,6,5,5), lwd= 1)
plot(mtcars$mpg, mtcars$hp,
     col = colours[gear_column],
     pch = 16,
     cex = 2, 
     xlab = "Miles per gallon",
     ylab = "Horsepower",
     cex.lab = 1)
legend(legend = sort(unique(mtcars$gear)),
       x = "right",
       fill = c("red", "green", "blue"))

#Plot the two plots that you have just made side-by-side by changing the global graphical parameters


#------------------

#in coding regions files, add column names, and add a new column containing the length of each region (you should have done this in the base R practical) 

#import file 
coding_gene_region <- read.table("data/coding_gene_region.bed") 

getwd()
colnames(coding_gene_region) <- c("Chromosome", "Start", "End", "ID", "Score", "Strand")
coding_gene_region$Genomic_interval_length <- coding_gene_region$End - coding_gene_region$Start

#Plot a histogram of the lengths using ggplot2: 

ggplot(coding_gene_region, aes(x = Genomic_interval_length)) +
        geom_histogram(bins = 100, colour = "red", fill = "blue") +
        labs(title = "Histogram of genomic intervals",
             x = "Genomic Interval Length",
             y = "Count") +
        theme(axis.title.y = element_text(size = 10, face = "bold"),
              axis.title.x = element_text(size = 20, face = "bold"),
              axis.text.x =  element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(angle = 45, hjust = 1, vjust = 0.5),
              plot.title = element_text(hjust = 0.5, size = 10)) +
        xlim(0, 500000) +
        scale_x_continuous(labels=comma)



#tidyverse count file 
sample_table <- read_tsv("data/tidyverse/obds_sampletable.tsv")
count_table <- read_tsv("data/tidyverse/obds_countstable.tsv.gz")
View(sample_table)
View(count_table)

### `pivot_longer()` - convert to long format

processed_count_table <- count_table %>% 
        pivot_longer(-Geneid, names_to = "sample", values_to = "count")

#join with gene information to get mgi_symbol 
listMarts()
ensembl <- useMart("ensembl") #connected to a database with datasets, we want to list datasets avaialble so we use list datasets 
datasets <- listDatasets(ensembl)
head(datasets) 
View(datasets)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
View(ensembl)        

filters <- listFilters(ensembl)
filters[1:5,]
View(filters)

attributes = listAttributes(ensembl)
attributes[1:5,]
View(attributes)

matching <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
      values = unique(processed_count_table$Geneid),
      mart = ensembl) 
View(matching)

#joining tables to match geneid note 
processed_count_table <- count_table %>%
        left_join(matching, by = c("Geneid", "ensembl_gene_id"))

#seperate out sample_title column into separate columns as it had 4 sets of data in 1 column 
processed_sample_table <- sample_table %>% 
        separate(sample_title, c("gene", "genotype", "cell_type", "replicates"), sep = "_") %>%
        unite("genotype", gene, genotype, sep = "_") %>%
        dplyr::select(-c(read_count, library_layout))
View(processed_sample_table)

processed_joined <- processed_count_table %>%
        left_join(processed_sample_table, by = c("sample" = "Sample_accession"))
#
calculated <- processed_joined %>%
        group_by(sample) %>%
        mutate(total_count_per_sample = sum(count)) %>%
        mutate(total_count_in_millions = total_count_per_sample / 1000000) %>%
        mutate(CPM = count / total_count_in_millions)  %>%
        mutate(logCPM = log2(CPM+1))
View(calculated)
