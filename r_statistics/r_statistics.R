#import libraries
library(ggplot2)
library(tidyverse)


#Generate a vector of 1000 normally distributed values with mean 10 and standard deviation 5.
num_vector <- rnorm(1000, mean = 10, sd = 5)

#Inspect the output of the summary()  function for that vector.
summary(num_vector)

#Compute the mean and standard deviation for those values

mean_vector <- mean(num_vector)
sd(num_vector)

#Compute the deciles (i.e. 10 evenly spaced quantiles) for those values.
quantile(num_vector, probs = seq(0, 1, 0.1))
#also can use: quantile(num_vector, probs = seq(0,1, length.out=10))

#Visualise the distribution of those values as a histogram. You can use base RR ggplot
#base R
hist(num_vector, breaks = 20)

#ggplot - need to convert to DF for ggplot
vector_df <- data.frame(NAME = num_vector)
head(vector_df)
ggplot(vector_df, aes(x= NAME)) +
    geom_histogram()

#Visualise as vertical lines on the histogram: the mean, median, one standard deviation, and one #median absolute deviation.'''

ggplot(vector_df, aes(x= NAME)) +
    geom_histogram() +
    geom_vline(xintercept = mean_vector)

#Generate a new vector with a lot more values (e.g., one million)
large_vector <-  rnorm(1e6, mean = 10, sd = 5)
summary(large_vector)
hist(num_vector, breaks = 100)

#more datapoint you have the more normally distributed it will look

#-----------------------------------#
#set seed as todays date
set.seed(20201014)

#make table with various distrubutions
dst <- tibble(rnorm=rnorm(1000))
#-----------------------------------#

#Plot the cumulative distribution in the range [-5,5]
x_index = seq(-5, 5, length.out = 1000)

cumulat_dist <- pnorm(q= seq(-5,5,0.01), mean = 0, sd = 1)
View(cumulat_dist)

dst <- tibble(x_coordinate = x_index, cdf_normal_distribution = pnorm(q = x_index, mean = 0, sd = 1))

dst %>%
    ggplot(aes(x=x_coordinate, y=cdf_normal_distribution)) +
    geom_point()

#Plot the inverse cumulative distribution function for quantiles in 0.01 increment
x_prob <- seq(0, 1, 0.01)

dist <- tibble(x_prob = x_prob, inverse_cdf_normal = qnorm(p = x_prob, mean = 0, sd = 1))

dist %>%
    ggplot(aes(x=x_prob, y=inverse_cdf_normal)) +
    geom_point()

#Plot the density function in the range [-5, 5]
dst <- tibble(x_coordinate = x_index, density_normal_distribution = dnorm(q = x_index, mean = 0, sd = 1))

dist %>%
    ggplot(aes(x=x_coorindate, y=density_normal_distribution)) +
    geom_point()


#What is the probability of observing a value greater than 2?
pnorm(2)
1 - pnorm(2)

#What is the probability of observing a value between -2 and 2?
1 - (pnorm(2) - pnorm(-2))



#ecdf
iris
ecdf(iris$Sepal.Length)

out <- ecdf(iris$Sepal.Length)
plot(out)
knots(out)


#---------------------------------


#Iris data set gives the measurements in centimeters of the variables sepal length and width and petal length and width, respectively, for 50

#1. Use the summary  function to view some information about each column.
summary(iris$Sepal.Length)
summary(iris$Species)
summary(iris$Sepal.Width)
summary(iris)

#2. Visualise the distribution of sepal length, stratified by species
ggplot(iris, aes(x= Sepal.Length), bins = 30) +
    geom_histogram(colour = "black") +
    facet_wrap(~Species, ncol = 1)

#3. Is sepal.length  length normally distributed? Overall? Within each species?

shapiro.test(iris$Sepal.Length)
#not normally distributed as p = 0.01018 which is larger than p=0.05 so we reject null hypothesis


species_setosa <- subset(iris, Species == "setosa", select = Sepal.Length)
species_setosa

species_versicolor <- subset(iris, Species == "versicolor", select = Sepal.Length)
species_versicolor

species_virginica <- subset(iris, Species == "virginica", select = Sepal.Length)
species_verginica

#or you could select by species and sepal length using group_by

#iris %>%
#group_by(Species) %>%
    #summarise(shapiro=shapiro.test(Sepal.Length)$p.value)

#or using a for loop
#for (species_name in as.character(unique(iris$Species))){​​
    #print(shapiro.test(iris$Sepal.Length[iris$Species==as.character(species_name)]))
    #}​​



#shaprio test requires vector so still need $<col> here
shapiro.test(species_setosa$Sepal.Length)
shapiro.test(species_versicolor$Sepal.Length)
shapiro.test(species_virginica$Sepal.Length)


#4. Is there a significant variation between sepal.length between various species

anova <- aov(formula = Sepal.Length ~ Species, data = iris)
summary(anova)

levels(iris$Species)
pair_species <- subset(iris, Species == "setosa", select = Sepal.Length)

ttest1 <- t.test(species_setosa$Sepal.Length, species_versicolor$Sepal.Length)
str(ttest1)

ttest2 <- t.test(species_setosa$Sepal.Length, species_virginica$Sepal.Length)
str(ttest2)

ttest3 <- t.test(species_virginica$Sepal.Length, species_versicolor$Sepal.Length)
str(ttest3)


#------------------#

#Fit a linear model to measure the effect of time & diet on weight
ChickWeight
head(ChickWeight)


#first formula is running 2 separate univariate linear model,

formula <- formula(weight ~ Diet + Time)
lm_chick <- lm(formula = formula, data = ChickWeight)
summary(lm_chick)
#coeffi estiamte is an average

ggplot(ChickWeight, aes(x= Time, y = weight, colour = Diet)) +
    geom_point() +
    geom_smooth(method = "lm") +
    geom_abline(slope = 8.7505, intercept = 10.9244)
#abline drawring average across diets

#taking diet and time instead of separately, now we're running multivariate
formula2 <- formula(weight ~ Diet * Time)
lm_chick <- lm(formula = formula2, data = ChickWeight)
summary(lm_chick)

#on avg all gain 6.8g over time. diet2 gives extra 1.7g on top etc etc

ChickWeight$Diet <- relevel(ChickWeight$Diet, "3")
lm_chick <- lm(formula = formula2, data = ChickWeight)
summary(lm_chick)
levels(ChickWeight$Diet) #all diets defined by diet 3 (reference)

ggplot(ChickWeight, aes(x= Time, y = weight, colour = Diet)) +
    geom_point() +
    geom_smooth(method = "lm") +
    geom_abline(slope = 11.4229, intercept = 18.2503) +
    geom_abline(slope = (11.4229 - 4.5811), intercept = (18.2503 + 12.6807))




#-------------------#

#setting row.names = says the first row is the header?
log_counts <- read.csv("data/logcounts.csv", row.names = 1)
View(log_counts)
cell_md <- read.csv("data/cell_metadata.csv",row.names=1)
View(cell_md)

#For each gene (i.e. row) in log_count, use cell_md and a statistical test of your choice to identify gene differentially expressed in cells infected with Salmonella relative to the control uninfected cells.

#
gene1 <- data_frame("log_count" = as.numeric(log_counts[1,]), infection = cell_md$Infection)
test_result <- t.test(log_count~infection, gene1)
str(test_result)
test_result[['p.value']]
#
diff_exp <- function(gene_index, matrix, groups){
    gene_row <- data.frame("log_count" = as.numeric(matrix[gene_index,]), infection = groups)
    test_result<- t.test(log_count~infection, gene_row)
    return(test_result[['p.value']])
}

diff_exp(3,log_counts,cell_md$Infection)

p_values <- vapply(seq(1,nrow(log_counts)), diff_exp, numeric(1), matrix = log_counts, groups = cell_md$Infection)

names(p_values) <- rownames(log_counts)
head(p_values)
p_adjusted <- p.adjust(p_values, method = "BH") #p-values adjusted with Benjamini & Hochberg - fewer p values
hist(p_adjusted)

head(sort(p_adjusted))

hist(p_values, breaks = 50)
hist(p_adjusted, breaks = 50)

table(p_adjusted < 0.5) #tells us how many true and false

#subset the significant genes and use them to find the signalling pathway
#which adjusted p values are less than 0.5
#give names of genes with p-value less than 0.5
sig_gene_names <- names(p_adjusted[p_adjusted < 0.05]) #extract out names of p_adjusted values less than 0.5
sig_gene_names




#--------------------#
#Fishers exact test:
#Fisher test - for each pathway we take list of genes (group 2) compared with group 1 - genes that are differntially expressed . Total genes = entire set of genes we've tested (1000 in this dataset but in real case scenario it's all the genes that could be differentially expressed - rarely use entire genome )

#group 1 = deferentially expressed gene
#group 2 = genes in your pathway

#Fishers more precise than Chi-sqr with small dataset -don't run chi-srq when expected values < 5 - always run Fishers

#build contingency table:

#read in human_go
human_go <- read.csv("data/human_go_bp.csv")
View(human_go)

#put as list
pathway_list <- split(human_go$ensembl_gene_id, human_go$go_id) #create list so column ensemble id is split and all genes will belong to 1 go id

path_1 <- pathway_list[[1]] #[[]] returns item in list, [] returns list

#build a table that has 1 row per 1000 genes (table with 1000 rows), each row represents a gene, 1st col (true/false - ), 2nd col (true,false)

sum(rownames(log_counts) %in% sig_gene_names) #for each gene in rownames(log_count) is it in significant_gene_names? -returns boolean

overlap_table <- data.frame(row.names = rownames(log_counts),
           significant = rownames(log_counts) %in% sig_gene_names,
           pathway = rownames(log_counts) %in% path_1)

contigency_table <- table(overlap_table) #used a dataframe to get the contigency table

#Fisher table:
fisher_result <- fisher.test(contigency_table, alternative = "less")
fisher_result$p.value

#fisher asks if pathway over or under represented in list? -over-representation uses alternative = "greater" in fisher.test gives p-value =1, under-representation uses alternative = "less"



#For all pathways:

#define a function for fisher.test, in order to run the signal pathway
signal_pathway <- function(pathway, matrix, sig_genes){
    overlap_table <- data.frame(
        row.names = rownames(matrix),
        significant = factor(rownames(matrix) %in% sig_genes, c(FALSE, TRUE)),
        pathway = factor(rownames(matrix) %in% pathway, c(FALSE, TRUE)))
    fisher_result <- fisher.test(table(overlap_table), alternative = "greater")
    return(fisher_result$p.value)
}


#significant = factor(rownames(matrix) %in% sig_genes, c(FALSE, TRUE)) - defines level of factor, even if R only sees 1 value, we're forcing R to create two levels (rows) to maintain contingency table structure.

#calling function on first gene pathway, defined earlier as path1
signal_pathway(path_1, log_counts, sig_gene_names)

#call function on all pathways using lapply

fisher_all_pathway <- sapply(pathway_list, signal_pathway, matrix = log_counts, sig_genes = sig_gene_names)

unique(fisher_all_pathway)

hist(fisher_all_pathway, breaks = 50)

fisher_final <-p.adjust(fisher_all_pathway, method = "BH")
hist(fisher_final, breaks = 50)
head(sort(fisher_final))


#bioconductor package for gene ontology analysis - topGO - plots graph with entire branches that are significant (overlapping pathways) - look for more generic pathway with lowest p-value. Report significant enrichment in the most generic pathway (don't list the 10 associated pathways just say for eg T cell proliferation)




#vapply - for each pathway - build table for features test (diff or not diff expressed) - how many genes belong in pathway or not
#Fishers test

#take first pathway, test 1 pathway - turn into function and test vapply over all pathways

#initialise env to new place
