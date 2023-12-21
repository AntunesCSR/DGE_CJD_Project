# AUTHOR: Baiba Vilne
# DATE:  2023-05-11

# DESCRIPTION:
# The code loads the limma and edgeR packages and defines the main() function. 
# The function reads in input files, combines the counts for the samples and 
# groups them by treatment, creates a design matrix for the linear model, 
# performs TMM normalization, calculates log-transformed CPM values, fits 
# a linear model, assigns gene symbols to the fitted model object, defines 
# the contrast matrix and fits the linear model with contrasts, calculates 
# moderated t-statistics using an empirical Bayes approach, extracts the top 
# differentially expressed genes, and writes them to a file. 

# Load required packages
library(limma)
library(edgeR)
library(dplyr)

# Read in the input files
countTable <- read.delim('gene_expr_data.csv', stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
#Remove outlier sample
countTable <- select(countTable, -CJD_CB_20)

#Extract clinical information to assign groups to samples
snames <- colnames(countTable)
split_names <- strsplit(snames, "_")
sample_group <- sapply(split_names, function(x) x[1])
sample_tissue <- sapply(split_names, function(x) x[2])
groups <- interaction(sample_group, sample_tissue)
groups

# Create a design matrix for the linear model
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
fit <- lmFit(countTable, design)

# Define the contrast matrix and fit the linear model with contrasts
contrast_matrix <- makeContrasts(
  CJD_vs_CT = (CJD.CB + CJD.FC) - (CT.CB + CT.FC),
  levels = design
)

# Fit the linear model
fit2 <- contrasts.fit(fit, contrast_matrix)

# Calculate moderated t-statistics using an empirical Bayes approach
eBayesFit <- eBayes(fit2, trend = TRUE)

# Extract the top differentially expressed genes
allGenes <- topTable(eBayesFit, number = Inf)
print(head(allGenes))

# Filter based on adjusted p-value and fold change
allGenes <- allGenes %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)

write.table(allGenes, file = "DEG_controls_vs_disease_filtered.csv", sep = "\t", quote = FALSE, row.names = TRUE)

###limma multiple groups

#Defining contrast matrix for analysis (4 groups)
contrast_matrix_ANOVA <- makeContrasts(
  CJD.CB_vs_CT.CB = CJD.CB - CT.CB,
  CJD.FC_vs_CT.FC = CJD.FC - CT.FC,
  CJD.CB_vs_CJD.FC = CJD.FC - CJD.CB,
  levels = design
)

#Performing limma
fit2_ANOVA <- contrasts.fit(fit, contrast_matrix_ANOVA)
eBayesFit_ANOVA <- eBayes(fit2_ANOVA)
top_table <- topTable(eBayesFit_ANOVA, number=Inf)
#Filter based on adjusted p-value
top_table <- top_table %>%
  filter(adj.P.Val < 0.05)

#Create table for cerebellum region
CB_DEG <- select(top_table, -CJD.FC_vs_CT.FC, -CJD.CB_vs_CJD.FC)
CB_DEG <- CB_DEG %>%
  filter(abs(CJD.CB_vs_CT.CB) > 1)
write.table(CB_DEG, file = "DEG_cerebellum.csv", sep = "\t", quote = FALSE, row.names = TRUE)

# Create table for frontal cortex region
FC_DEG <- select(top_table, -CJD.CB_vs_CT.CB, -CJD.CB_vs_CJD.FC)
FC_DEG <- FC_DEG %>%
  filter(abs(CJD.FC_vs_CT.FC) > 1)
write.table(FC_DEG, file = "DEG_cortex.csv", sep = "\t", quote = FALSE, row.names = TRUE)

# Create table comparing region differences
regions_DEG <- select(top_table, -CJD.CB_vs_CT.CB, -CJD.FC_vs_CT.FC)
regions_DEG <- regions_DEG %>%
  filter(abs(CJD.CB_vs_CJD.FC) > 1)
write.table(regions_DEG, file = "DEG_regions.csv", sep = "\t", quote = FALSE, row.names = TRUE)
