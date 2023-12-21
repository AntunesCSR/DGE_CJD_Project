library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dplyr)

#Data table preparation: loading necessary data, excluding outlier and renaming row for clarity
DEG_CB <- read.delim('DEG_cerebellum.csv', stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
DEG_FC <- read.delim('DEG_cortex.csv', stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
DEG_dis <- read.delim('DEG_regions.csv', stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)

colnames(DEG_CB)[colnames(DEG_CB) == "row.names"] <- "genes"
colnames(DEG_FC)[colnames(DEG_FC) == "row.names"] <- "genes"
colnames(DEG_dis)[colnames(DEG_dis) == "row.names"] <- "genes"

genes_expressed <- read.delim('gene_expr_data.csv', stringsAsFactors = FALSE, check.names = FALSE, header = TRUE, row.names = NULL)
colnames(genes_expressed)[colnames(genes_expressed) == "row.names"] <- "gene_symbol"
genes_expressed <- select(genes_expressed, -CJD_CB_20)

#Selecting only relevant samples for each heatmap
snames <- colnames(genes_expressed)
snames_1 <- snames[-1]

selected_cb_columns <- snames_1[grep("CB", snames_1)]
selected_fc_columns <- snames_1[grep("FC", snames_1)]
selected_dis_columns <- snames_1[grep("CJD", snames_1)]

selected_columns_CB <- c(snames[1], selected_cb_columns)
selected_columns_FC<- c(snames[1], selected_fc_columns)
selected_columns_dis <- c(snames[1], selected_dis_columns)

selected_data_CB <- genes_expressed[, selected_columns_CB]
selected_data_FC <- genes_expressed[, selected_columns_FC]
selected_data_dis <- genes_expressed[, selected_columns_dis]

#Selecting differentially expressed genes from the gene expression matrix
i_CB <- which(selected_data_CB$gene_symbol %in% DEG_CB$genes)
i_FC <- which(selected_data_FC$gene_symbol %in% DEG_FC$genes)
i_dis <- which(selected_data_dis$gene_symbol %in% DEG_dis$genes)

#Transferring the data table into matrix
numeric_data_CB <- selected_data_CB[, sapply(selected_data_CB, is.numeric)]
numeric_data_FC <- selected_data_FC[, sapply(selected_data_FC, is.numeric)]
numeric_data_dis <- selected_data_dis[, sapply(selected_data_dis, is.numeric)]

genes_expressed[, sapply(genes_expressed, is.numeric)] <- lapply(genes_expressed[, sapply(genes_expressed, is.numeric)], as.numeric)

genes_expressed_matrix_CB <- data.matrix(numeric_data_CB)
genes_expressed_matrix_FC <- data.matrix(numeric_data_FC)
genes_expressed_matrix_dis <- data.matrix(numeric_data_dis)

#Selecting column names for heatmap labels
column_names_CB <- colnames(numeric_data_CB)
column_names_FC <- colnames(numeric_data_FC)
column_names_dis <- colnames(numeric_data_dis)

#Z-score normalisation
gene_expr_matrix_CB <- t(scale(t(genes_expressed_matrix_CB)))
gene_expr_matrix_FC <- t(scale(t(genes_expressed_matrix_FC)))
gene_expr_matrix_dis <- t(scale(t(genes_expressed_matrix_dis)))

#Extracting groups for heatmap annotation
snames_2_CB <- colnames(numeric_data_CB)
snames_2_FC <- colnames(numeric_data_FC)
snames_2_dis <- colnames(numeric_data_dis)

split_names_CB <- strsplit(snames_2_CB, "_")
split_names_FC <- strsplit(snames_2_FC, "_")
split_names_dis <- strsplit(snames_2_dis, "_")

sample_group_CB <- sapply(split_names_CB, function(x) x[1])
sample_group_FC <- sapply(split_names_FC, function(x) x[1])
sample_group_dis <- sapply(split_names_dis, function(x) x[2])

#Contruct heatmap annotations based on extracted groups. TODO: check why colours are not working...
ha_CB = HeatmapAnnotation(
  Group = sample_group_CB, 
  col = list(Group = c("CJD"="#e76f51", "CT"="#287271")
  )
)
ha_FC = HeatmapAnnotation(
  Group = sample_group_FC, 
  col = list(Group = c("CJD"="#e76f51", "CT"="#287271")
  )
)
ha_dis = HeatmapAnnotation(
  Group = sample_group_dis, 
  col = list(Group = c("FC" = "#ff9f1c", "CB" = "#ffbf69")
  )
)

#Plotting the heatmaps. TODO: add name for the legend (Z-score)
png("heatmap_DEG_CB.png", width = 3000, height = 3000, res = 300)
Heatmap(gene_expr_matrix_CB[i_CB,], 
        show_row_names = FALSE,
        cluster_rows = TRUE,
        top_annotation = ha_CB,
        name = "Z-score") 
dev.off()
png("heatmap_DEG_FC.png", width = 3000, height = 3000, res = 300)
Heatmap(gene_expr_matrix_FC[i_FC,], 
        show_row_names = FALSE,
        cluster_rows = TRUE,
        top_annotation = ha_FC,
        name = "Z-score") 
dev.off()
png("heatmap_DEG_REG.png", width = 3000, height = 3000, res = 300)
Heatmap(gene_expr_matrix_dis[i_dis,], 
        show_row_names = FALSE,
        cluster_rows = TRUE,
        top_annotation = ha_dis,
        name = "Z-score") 
dev.off()


