# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(readxl)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x <- stringr::str_replace(stringr::str_replace(stringr::str_replace(x,"MRNA", "mRNA"), "NcRNA", "ncRNA"), "RRNA", "rRNA")
  return(x)
}


# Set working directory
setwd("PATH/TO/DIRECTORY")

# Load top tables
load("Output/topTable.RData")

# Add annotation to top table
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
annotations <- getBM(attributes=c("ensembl_gene_id",
                                  "hgnc_symbol"), 
                     filters = 'ensembl_gene_id',
                     values = topTable_complex$ID,
                     mart = ensembl)

# Add gene symbol to top table
topTable_complex <- left_join(topTable_complex, annotations, by = c("ID" = "ensembl_gene_id"))
topTable_KO <- left_join(topTable_KO, annotations, by = c("ID" = "ensembl_gene_id"))
topTable_WT <- left_join(topTable_WT, annotations, by = c("ID" = "ensembl_gene_id"))


# Get downregulated genes
down_complex <- topTable_complex$ID[(topTable_complex$logFC < 0) & (topTable_complex$FDR < 0.05)]
down_WT <- topTable_WT$ID[(topTable_WT$logFC < 0) & (topTable_WT$FDR < 0.05)]
down_KO <- topTable_KO$ID[(topTable_KO$logFC < 0) & (topTable_KO$FDR < 0.05)]

up_complex <- topTable_complex$ID[(topTable_complex$logFC > 0) & (topTable_complex$FDR < 0.05)]
up_WT <- topTable_WT$ID[(topTable_WT$logFC > 0) & (topTable_WT$FDR < 0.05)]
up_KO <- topTable_KO$ID[(topTable_KO$logFC > 0) & (topTable_KO$FDR < 0.05)]



#==============================================================================#
# Complex comparison
#==============================================================================#

geneList <- sort(setNames(topTable_complex$logFC,
                          topTable_complex$ID), decreasing = TRUE)

GOresults <- gseGO(geneList     = geneList,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   keyType = "ENSEMBL",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = Inf,
                   verbose      = FALSE)

GOtable_complex <- GOresults@result


#==============================================================================#
# WT comparison
#==============================================================================#

geneList <- sort(setNames(topTable_WT$logFC,
                          topTable_WT$ID), decreasing = TRUE)

GOresults <- gseGO(geneList     = geneList,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   keyType = "ENSEMBL",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = Inf,
                   verbose      = FALSE)

GOtable_WT <- GOresults@result


#==============================================================================#
# KO comparison
#==============================================================================#

geneList <- sort(setNames(topTable_KO$logFC,
                          topTable_KO$ID), decreasing = TRUE)

GOresults <- gseGO(geneList     = geneList,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   keyType = "ENSEMBL",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = Inf,
                   verbose      = FALSE)

GOtable_KO <- GOresults@result


save(GOtable_complex, GOtable_KO, GOtable_WT, file = "Output/GOtable.RData")


clusterProfiler::simplify


# Make plot
load("Output/GOtable.RData")

# Make similarity matrix of GO terms
simMatrix <- calculateSimMatrix(
  GOtable_complex$ID[1:100],
  orgdb="org.Hs.eg.db",
  ont = "BP",
  method = "Resnik"
)
save(simMatrix, file = "Output/simMatrix.RData")


scores <- setNames(-log10(GOtable_complex$pvalue)[1:100], GOtable_complex$ID[1:100])
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")



parent <- unique(reducedTerms$parent)
GOtable_red <- GOtable_complex[GOtable_complex$ID %in% parent,]

selGOterms <- GOtable_red$ID[1:10]

plotDF <- data.frame(GO = selGOterms,
                     Description = factor(rep(firstup(GOtable_complex[selGOterms, "Description"]),2),
                                          levels = rev(firstup(GOtable_complex[selGOterms, "Description"]))),
                     NES = c(GOtable_WT[selGOterms, "NES"],
                             GOtable_KO[selGOterms, "NES"]),
                     logP = c(-log10(GOtable_WT[selGOterms, "pvalue"]),
                              -log10(GOtable_KO[selGOterms, "pvalue"])),
                     Group = rep(c("WT", "KO"), each = length(selGOterms)))



GOplot <- ggplot(plotDF[plotDF$Group == "WT",]) +
  geom_bar(aes(x = logP, y = Description, fill = NES),
           position = position_dodge(), stat = "identity") +
  scale_fill_gradient2(low = "#2171B5", 
                       mid = "white", 
                       high = "#CB181D", 
                       midpoint = 0,
                       limits = c(-2.5,2.5),
                       oob = scales::squish) +
  xlab(expression(-log[10]~"p-value")) +
  ylab(NULL) +
  xlim(0,10) +
  ggtitle("PAR1 WT: S vs C") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(GOplot, width = 6, height = 4, file = "Output/GOplot_WT.png")


GOplot <- ggplot(plotDF[plotDF$Group == "KO",]) +
  geom_bar(aes(x = logP, y = Description, fill = NES),
           position = position_dodge(), stat = "identity") +
  scale_fill_gradient2(low = "#2171B5", 
                       mid = "white", 
                       high = "#CB181D", 
                       midpoint = 0,
                       limits = c(-2.5,2.5),
                       oob = scales::squish) +
  xlab(expression(-log[10]~"p-value")) +
  xlim(0,10) +
  ylab(NULL) +
  ggtitle("PAR1 KO: S vs C") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(GOplot, width = 6, height = 4, file = "Output/GOplot_KO.png")





