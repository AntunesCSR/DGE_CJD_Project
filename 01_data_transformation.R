#Transform the data into readable matrix for further analysis
dat <- read.delim("series_matrix_RTU.txt", stringsAsFactors=F, check.names=F)
counts <- as.matrix(dat[, -1])
rownames(counts) <- dat$ID_REF

# Load gene expression matrix 
gene_expression_matrix <- counts

# Load clinical information
clinical_info_table <- read.table("sample_data.txt", header = TRUE)

# Create a mapping of sample names to new column names
new_column_names <- apply(clinical_info_table, 1, function(row) {
  paste(row["Diagnosis"], row["SampleType"], row["PatientNumber"], sep = "_")
})

# Rename the columns in the gene expression matrix
colnames(gene_expression_matrix) <- new_column_names

# Print the updated matrix with new column names
head(gene_expression_matrix)

write.table(gene_expression_matrix, file = "gene_expr_data.csv", sep = "\t", quote = FALSE, row.names = TRUE)

