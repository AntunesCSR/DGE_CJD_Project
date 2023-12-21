library(ggfortify)
library(ggbiplot)
library(devtools)
library("factoextra")

#Extract the data and group names for PCA plot
dat <- read.delim("gene_expr_data.csv", stringsAsFactors=F, check.names=F)
snames <- colnames(dat)
split_names <- strsplit(snames, "_")
sample_group <- sapply(split_names, function(x) x[1])
sample_tissue <- sapply(split_names, function(x) x[2])
groups <- interaction(sample_group, sample_tissue)
groups

#Transform into matrix for PCA 
dat <- t(dat)

#Calculate PCA and create the plot
dat.pca <- prcomp(dat, center = TRUE, scale = FALSE)
summary(dat.pca)

# Plot PCA - Original
# jpeg("rplot_PCA.jpg", width = 350, height = 350)
# fviz_pca_ind(dat.pca,
#              geom.ind = "point",
#              col.ind = groups, 
#              addEllipses = TRUE,
#              ellipse.level=0.85,
#              legend.title = "Groups"
#              ) +
#              labs(title = "PCA plot", x = "PC1 (47.8%)", y = "PC2 (24.2%)")
# dev.off()

# Plot PCA - Modification 
# Set up a color palette for better visualization
my_palette <- c("#e76f51", "#287271", "#f4a261", "#2a9d8f", "#984EA3", "#FFC20A")

# Increase the size of points and use a better color palette
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point",
                         col.ind = groups,
                         palette = my_palette,
                         pointsize = 3,  # Adjust the size of points
                         addEllipses = TRUE,
                         ellipse.level = 0.85,
                         legend.title = "Groups"
) +
  theme_light() + # Light theme for better visualization
  theme(legend.position = "right") +  # Adjust the legend position
  labs(title = "Principal Component Analysis (PCA) Plot",
       subtitle = "Expression Data",
       #caption = "Source: Your Data Source",
       x = "PC1 (47.8%)",
       y = "PC2 (24.2%)"
  ) +
  theme(plot.title = element_text(size = 18, face = "bold"),  # Adjust title font size and style
        plot.subtitle = element_text(size = 14, color = "gray"),  # Adjust subtitle font size and color
        plot.caption = element_text(size = 12, color = "darkgray"),  # Adjust caption font size and color
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")
  )

# Save the plot as a PNG file
ggsave("rplot_PCA.png", pca_plot, width = 10, height = 10, units = "in", dpi = 300)

