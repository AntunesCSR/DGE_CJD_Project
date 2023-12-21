library("VennDiagram")  
library("ggvenn") 

data.CB <- read.delim('DEG_cerebellum.csv', stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
colnames(data.CB)[colnames(data.CB) == "row.names"] <- "genes"
data.FC <- read.delim('DEG_cortex.csv', stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
colnames(data.FC)[colnames(data.FC) == "row.names"] <- "genes"
CB.genes <- data.CB$genes
FC.genes <- data.FC$genes

gene_lists <- list('CB-exclusive sCJD DEGs' = CB.genes, 'FC-exclusive sCJD DEGs' = FC.genes)
ggvenn(gene_lists)
#TODO: change colours and fonts

