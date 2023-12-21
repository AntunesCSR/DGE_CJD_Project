######################### Installations required to run the scripts in this repository #########################

# Install and load BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of Bioconductor packages
bioconductor_packages <- c("edgeR", "genefilter", "limma", "pheatmap", "DESeq2")

# Install Bioconductor packages
for (package in bioconductor_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

# List of CRAN packages
cran_packages <- c("ggfortify", "ggbiplot", "devtools", "R.utils", "pbapply")

# Install CRAN packages
for (package in cran_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

