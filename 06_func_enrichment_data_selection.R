library("dplyr")

func_KEGG <- read.delim("chart_KEGG.txt")
func_KEGG <- func_KEGG[1:40,]
write.table(func_KEGG, file = "func_enrich_KEGG.csv", sep = "\t", quote = FALSE, row.names = TRUE)

func_GO <- read.delim("chart_GO.txt")
func_GO <- func_GO[1:40,]
write.table(func_GO, file = "func_enrich_GO.csv", sep = "\t", quote = FALSE, row.names = TRUE)

func_BIOCARTA <- read.delim("chart_BIOCARTA.txt")
func_BIOCARTA <- func_BIOCARTA[1:40,]
write.table(func_BIOCARTA, file = "func_enrich_BIOCARTA.csv", sep = "\t", quote = FALSE, row.names = TRUE)

