# Functional Analysis Term plots 

# Example for loading text files
chart_BIOCARTA <- read.table("chart_BIOCARTA.txt", header = TRUE, sep = "\t")
chart_GO <- read.table("chart_GO.txt", header = TRUE, sep = "\t")
chart_KEGG <- read.table("chart_KEGG.txt", header = TRUE, sep = "\t")

# Example for loading CSV files
func_enrich_BIOCARTA <- read.csv("func_enrich_BIOCARTA.csv", header = TRUE, sep = "\t")
func_enrich_GO <- read.csv("func_enrich_GO.csv", header = TRUE, sep = "\t")
func_enrich_KEGG <- read.csv("func_enrich_KEGG.csv", header = TRUE, sep = "\t")

# Combine data frames for each database
combined_BioCarta <- rbind(chart_BIOCART, func_enrich_BIOCARTA)
combined_GO <- rbind(chart_GO, func_enrich_GO)
combined_KEGG <- rbind(chart_KEGG, func_enrich_KEGG)

# Check number of common elements
length(intersect(combined_BioCarta$Term, combined_GO$Term))
length(intersect(combined_BioCarta$Term, combined_KEGG$Term))
length(intersect(combined_GO$Term, combined_KEGG$Term))


# Set a p-value threshold
pvalue_threshold <- 0.05

# Filter data frames based on the p-value threshold
filtered_GO <- func_enrich_GO[func_enrich_GO$PValue < pvalue_threshold, ]
filtered_KEGG <- func_enrich_KEGG[func_enrich_KEGG$PValue < pvalue_threshold, ]
filtered_BioCarta <- func_enrich_BIOCARTA[func_enrich_BIOCARTA$PValue < pvalue_threshold, ]

# Extract relevant columns for the filtered data
pathways_GO <- filtered_GO$Term
gene_counts_GO <- filtered_GO$Count

pathways_KEGG <- filtered_KEGG$Term
gene_counts_KEGG <- filtered_KEGG$Count

pathways_BioCarta <- filtered_BioCarta$Term
gene_counts_BioCarta <- filtered_BioCarta$Count

#### GO PLOT ###
# Save the bar plot as a PNG file
png("barplot_GO.png", width = 3000, height = 3000, res = 300)
par(mar = c(20, 6, 6, 4) + 0.1)  # Increase the bottom margin for more space

# Create a bar plot with angled and staggered x-axis labels
bar_heights <- barplot(gene_counts_GO, las = 2, cex.names = 0.8, col = "#52796f", main = "Gene Ontology Pathways", ylab = "Gene Count", width = 1.0)  # Adjust width as needed

# Stagger the x-axis labels with increased rotation and adjustment
text(x = bar_heights, y = par("usr")[3] - 0.5, labels = pathways_GO, srt = 60, adj = c(1, 0.5), xpd = TRUE, cex = 0.7)

# Save the plot
dev.off()

#### KEGG PLOT ###
# Save the bar plot as a PNG file
png("barplot_KEGG.png", width = 3000, height = 3000, res = 300)
par(mar = c(15, 6, 6, 4) + 0.1)  # Increase the bottom margin for more space

# Create a bar plot with angled and staggered x-axis labels
bar_heights <- barplot(gene_counts_KEGG, las = 2, cex.names = 0.8, col = "#e76f51", main = "KEGG Pathways", ylab = "Gene Count", width = 1.0)  # Adjust width as needed

# Stagger the x-axis labels with increased rotation and adjustment
text(x = bar_heights, y = par("usr")[3] - 0.5, labels = pathways_KEGG, srt = 60, adj = c(1, 0.5), xpd = TRUE, cex = 0.7)

# Save the plot
dev.off()

#### BioCarta PLOT ###
# Save the bar plot as a PNG file
png("barplot_BioCarta.png", width = 3000, height = 3000, res = 300)
par(mar = c(25, 6, 6, 4) + 0.1)  # Increase the bottom margin for more space

# Create a bar plot with angled and staggered x-axis labels
bar_heights <- barplot(gene_counts_BioCarta, las = 2, cex.names = 0.8, col = "#ff9f1c", main = "BioCarta Pathways", ylab = "Gene Count", width = 1.0)  # Adjust width as needed

# Stagger the x-axis labels with increased rotation and adjustment
text(x = bar_heights, y = par("usr")[3] - 0.5, labels = pathways_BioCarta, srt = 60, adj = c(1, 0.5), xpd = TRUE, cex = 0.7)

# Save the plot
dev.off()




