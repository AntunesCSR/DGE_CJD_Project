#Note to group: this does not work, our analysis has too little genes. Used other tool and saved this as an example...

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

#load table with list of differentially expressed genes
de_genes <- read.delim("DEG_controls_vs_disease_filtered.csv")

#load table with all genes detected to be expressed in the experiment
detected_genes <- read.delim('gene_expr_data.csv', stringsAsFactors = FALSE, check.names = FALSE, header = FALSE)
detected_genes <- detected_genes[-1,]

univgenes=unique(detected_genes[,1])

#perform over-representation analysis of GO BP terms

ego <- enrichGO(gene          = rownames(de_genes),
                universe      = univgenes,
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.2,
                qvalueCutoff  = 0.1,
                minGSSize = 10,
                maxGSSize = 500,
                readable      = TRUE)

ego@result$FoldEnrichment=parse_ratio(ego@result$GeneRatio)/parse_ratio(ego@result$BgRatio)


#to check the significant results in a table
egodf=as.data.frame(ego)

#to remove enriched terms with small number of associated genes
todrop=ego@result$ID[ego@result$Count<4]
ego=dropGO(ego,term=todrop)

egodf=as.data.frame(ego)

#to view results in a plot

dotplot(ego, x="FoldEnrichment",color="p.adjust",size="Count",showCategory=20,orderBy="Count")


#to export the results table

write.table(egodf,"egodf.txt",row.names=F)


#perform over-representation analysis of KEGG pathways

#this analysis requires using a different gene id
gene.df <- bitr(gene= rownames(de_genes), fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

univ.df <- bitr(gene= univgenes, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene          = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pAdjustMethod = "BH",
                 universe= univ.df$ENTREZID,
                 pvalueCutoff = 0.05)

kk@result$FoldEnrichment=parse_ratio(kk@result$GeneRatio)/parse_ratio(kk@result$BgRatio)

dotplot(kk, x="FoldEnrichment",color="p.adjust",size="Count",showCategory=20,orderBy="Count")


kkdf=as.data.frame(kk)

#view results on kegg map
browseKEGG(kk, 'hsa04919')
