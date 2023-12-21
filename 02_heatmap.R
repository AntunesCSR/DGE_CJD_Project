library(ComplexHeatmap)

data <- read.delim("gene_expr_data.csv", stringsAsFactors=F, check.names=F)

#Select most variable genes based on standard deviation
gene_variability <- apply(data, 1, sd)

# Rank genes by variability
ranked_genes <- order(gene_variability, decreasing = TRUE)

# Select a threshold
num_selected_genes <- 100
selected_genes <- ranked_genes[1:num_selected_genes]

# Extract the most variable genes
most_variable_genes <- data[selected_genes, ]

m <- as.matrix(most_variable_genes)
m <- t(scale(t(m)))


#Extract group names for the heatmap
snames <- colnames(data)
split_names <- strsplit(snames, "_")
sample_group <- sapply(split_names, function(x) x[1])
sample_tissue <- sapply(split_names, function(x) x[2])
sample_tissue


#Construct heatmap annotations based on extracted groups - modify colors
ha = HeatmapAnnotation(
  Group = sample_group,
  Tissue = sample_tissue,
  col = list(
    Group = c("CJD" = "#2ec4b6", "CT" = "#cbf3f0"),
    Tissue = c("FC" = "#ff9f1c", "CB" = "#ffbf69")
  )
)

# #Create heatmap - without clustering
# png("rplot_heatmap.png", width = 3000, height = 3000, res=300)
# Heatmap(m, show_row_names = FALSE,
#         cluster_rows = FALSE,
#         top_annotation = ha,
#         name = " ")
# dev.off()


#Create heatmap - with clustering
png("rplot_heatmap_HC.png", width = 3000, height = 3000, res=300)
Heatmap(m, show_row_names = FALSE,
        cluster_rows = TRUE,
        top_annotation = ha,
        name = "Z-Score")
dev.off()



