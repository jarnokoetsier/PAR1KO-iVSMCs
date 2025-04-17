# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(biomaRt)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)

# Set working directory
setwd("PATH/TO/DIRECTORY")

# Load top tables
load("topTable.RData")

# Add annotation to top table
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
annotations <- getBM(attributes=c("ensembl_gene_id",
                                  "entrezgene_id",
                                  "hgnc_symbol"), 
                     filters = 'ensembl_gene_id',
                     values = topTable_complex$ID,
                     mart = ensembl)

# Select significant genes (FDR-threshold = 0.01)
sigGenes <- annotations$hgnc_symbol[annotations$ensembl_gene_id %in%
                                      topTable_complex$ID[topTable_complex$FDR < 0.01]]

# Load GO-gene annotation
load("Data/GOgenes_BP_SYMBOL_Hs.RData")

# Select genes of interest:
selGenes <- GOgenes[["GO:0007169"]] # Cell surface receptor protein tyrosine kinase signaling pathway
selGenes <- GOgenes[["GO:0034599"]] # Response to oxidative stress

# Selected genes: in GO terms + statistically significant
selGenes <- intersect(selGenes, sigGenes)

# Prepare data for plotting
plotDF <- data.frame(ID = c(topTable_KO$ID, topTable_WT$ID),
                     logFC = c(topTable_KO$logFC, topTable_WT$logFC),
                     Group = c(rep("PAR1 KO: S vs C", nrow(topTable_KO)),
                               rep("PAR1 WT: S vs C", nrow(topTable_WT))))

# Add gene annotations to data frame
plotDF <- inner_join(plotDF, annotations, by = c("ID"= "ensembl_gene_id"))
plotDF <- unique(plotDF[,c("hgnc_symbol", "ID", "logFC", "Group")])
colnames(plotDF) <- c("Gene", "ID","logFC", "Group")

# Filter for genes of interest
plotDF_fil <- plotDF[plotDF$Gene %in% selGenes,]

# Perform clustering to get the order of the genes for the heatmap:

# 1. Make matrix
m <- as.data.frame(inner_join(topTable_KO[topTable_KO$ID %in% plotDF_fil$ID,
                                          c(6,1)], 
                              topTable_WT[topTable_WT$ID %in% plotDF_fil$ID,
                                          c(6,1)], by = "ID"))
rownames(m) <- m$ID
m <- as.matrix(m[,-1])

# 2. Set maximum absolute value to 3
m[m > 3]<- 3
m[m < -3] <- -3

# 3. Perform clustering (Euclidean distance and Ward D2 linkage methods)
model_genes <- hclust(dist(m,
                           method = "euclidean"), "ward.D2")

# Set order of genes
id_ordered <- model_genes$labels[model_genes$order]
plotDF_fil$ID <- factor(plotDF_fil$ID , levels = id_ordered)
genes_ordered <- arrange(plotDF_fil , by = ID)$Gene
plotDF_fil$Gene <- factor(plotDF_fil$Gene, levels = unique(genes_ordered))

# Make plot
p <- ggplot(plotDF_fil) +
  geom_tile(aes(y = Gene, x = Group, fill = logFC)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Transmembrane receptor protein\ntyrosine kinase signalling") +
  scale_fill_gradient2(low = "#2171B5", 
                       mid = "white", 
                       high = "#CB181D", 
                       midpoint = 0,
                       limits = c(-3,3),
                       oob = scales::squish) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# Save plot
ggsave(p, file = "tyrosine kinase_heatmap.png", width = 4, height = 12)

