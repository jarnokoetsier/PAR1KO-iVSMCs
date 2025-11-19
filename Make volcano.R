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
library(ggVennDiagram)

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


#==============================================================================#
# Complex comparison
#==============================================================================#

# Make volcano plot
topTable_sel <- topTable_complex
p = "raw"
p_threshold = 0.05
logFC_threshold = 0
unchanged_color = "darkgrey"
down_color = "#4292C6"
up_color = "#CB181D"

plotDF <- data.frame(log2FC = topTable_sel[,"logFC"],
                     Pvalue = topTable_sel[,"PValue"],
                     logP = -log10(topTable_sel[,"PValue"]),
                     adjPvalue = topTable_sel[,"FDR"],
                     logAdjP = -log10(topTable_sel[,"FDR"]),
                     GeneID = topTable_sel[,"ID"],
                     Name = topTable_sel[,"hgnc_symbol"])

plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$adjPvalue > p_threshold] <- "unchanged"
plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$adjPvalue <= p_threshold] <- "downregulated"
plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$adjPvalue <= p_threshold] <- "upregulated"


selGenes <- c("PRSS35", "PLAT", "DACT1", "MEDAG", "CDKN2B", "TFPI",
              "SPOCD1", "CYP1B1", "OLFM2", "COL21A1")
# selGenes <- c(topTable_sel$hgnc_symbol[topTable_sel$logFC < 0][1:5],
#               topTable_sel$hgnc_symbol[topTable_sel$logFC > 0][1:5])

set.seed(123)
volcano <- ggplot2::ggplot() +
  ggplot2:: geom_point(data = plotDF[plotDF$Colour == "unchanged",],
                       ggplot2::aes(x = log2FC, y = logP, color = "Unchanged",
                                    alpha = logP, size = logP)) +
  ggplot2::geom_point(data = plotDF[plotDF$Colour == "downregulated",],
                      ggplot2::aes(x = log2FC, y = logP, color = "Downregulated",
                                   alpha = logP, size = logP)) +
  ggplot2::geom_point(data = plotDF[plotDF$Colour == "upregulated",],
                      ggplot2::aes(x = log2FC, y = logP, color = "Upregulated",
                                   alpha = logP, size = logP),) +
  #ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "lightgrey") +
  ggrepel::geom_text_repel(data = plotDF[plotDF$Name %in% selGenes,],
                           ggplot2::aes(x = log2FC, y = logP,
                                        label = Name),
                           min.segment.length = 0.1, force_pull = 0.2) +
  ggplot2::scale_color_manual(values = setNames(c(unchanged_color, down_color, up_color),
                                                c("Unchanged", "Downregulated", "Upregulated"))) +
  ggplot2::xlab(expression(log[2]~"FC")) +
  ggplot2::ylab(expression(-log[10]~"P value")) +
  ggplot2::xlim(-10,10) +
  ggplot2::ylim(0,23) +
  ggtitle("(PAR1 KO: S vs C) vs (PAR WT: S vs C)") +
  guides(size = "none", alpha = "none") +
  ggplot2::labs(color = NULL) + 
  ggplot2::scale_size_continuous(range = c(1,3)) +
  ggplot2::scale_alpha_continuous(range = c(0.5,1)) +
  ggplot2::theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

volcano

ggsave(volcano, width = 6, height = 5, file = "Output/Volcano_complex.png")

#==============================================================================#
# WT comparison
#==============================================================================#

# Make volcano plot
topTable_sel <- topTable_WT
p = "raw"
p_threshold = 0.05
logFC_threshold = 0
unchanged_color = "darkgrey"
down_color = "#4292C6"
up_color = "#CB181D"

plotDF <- data.frame(log2FC = topTable_sel[,"logFC"],
                     Pvalue = topTable_sel[,"PValue"],
                     logP = -log10(topTable_sel[,"PValue"]),
                     adjPvalue = topTable_sel[,"FDR"],
                     logAdjP = -log10(topTable_sel[,"FDR"]),
                     GeneID = topTable_sel[,"ID"],
                     Name = topTable_sel[,"hgnc_symbol"])

plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$adjPvalue > p_threshold] <- "unchanged"
plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$adjPvalue <= p_threshold] <- "downregulated"
plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$adjPvalue <= p_threshold] <- "upregulated"


selGenes <- c("PRSS35", "PLAT", "DACT1", "MEDAG", "CDKN2B", "TFPI",
              "SPOCD1", "CYP1B1", "OLFM2", "COL21A1")
# selGenes <- c(topTable_sel$hgnc_symbol[topTable_sel$logFC < 0][1:5],
#               topTable_sel$hgnc_symbol[topTable_sel$logFC > 0][1:5])

set.seed(123)
volcano <- ggplot2::ggplot() +
  ggplot2:: geom_point(data = plotDF[plotDF$Colour == "unchanged",],
                       ggplot2::aes(x = log2FC, y = logP, color = "Unchanged",
                                    alpha = logP, size = logP)) +
  ggplot2::geom_point(data = plotDF[plotDF$Colour == "downregulated",],
                      ggplot2::aes(x = log2FC, y = logP, color = "Downregulated",
                                   alpha = logP, size = logP)) +
  ggplot2::geom_point(data = plotDF[plotDF$Colour == "upregulated",],
                      ggplot2::aes(x = log2FC, y = logP, color = "Upregulated",
                                   alpha = logP, size = logP),) +
  #ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "lightgrey") +
  ggrepel::geom_text_repel(data = plotDF[plotDF$Name %in% selGenes,],
                           ggplot2::aes(x = log2FC, y = logP,
                                        label = Name),
                           min.segment.length = 0.1, force_pull = 0.2) +
  ggplot2::scale_color_manual(values = setNames(c(unchanged_color, down_color, up_color),
                                                c("Unchanged", "Downregulated", "Upregulated"))) +
  ggplot2::xlab(expression(log[2]~"FC")) +
  ggplot2::ylab(expression(-log[10]~"P value")) +
  ggplot2::xlim(-10,10) +
  ggplot2::ylim(0,23) +
  ggtitle("PAR1 WT: S vs C") +
  guides(size = "none", alpha = "none") +
  ggplot2::labs(color = NULL) + 
  ggplot2::scale_size_continuous(range = c(1,3)) +
  ggplot2::scale_alpha_continuous(range = c(0.5,1)) +
  ggplot2::theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

volcano

ggsave(volcano, width = 6, height = 5, file = "Output/Volcano_WT_alt.png")


#==============================================================================#
# KO comparison
#==============================================================================#

# Make volcano plot
topTable_sel <- topTable_KO
p = "raw"
p_threshold = 0.05
logFC_threshold = 0
unchanged_color = "darkgrey"
down_color = "#4292C6"
up_color = "#CB181D"

plotDF <- data.frame(log2FC = topTable_sel[,"logFC"],
                     Pvalue = topTable_sel[,"PValue"],
                     logP = -log10(topTable_sel[,"PValue"]),
                     adjPvalue = topTable_sel[,"FDR"],
                     logAdjP = -log10(topTable_sel[,"FDR"]),
                     GeneID = topTable_sel[,"ID"],
                     Name = topTable_sel[,"hgnc_symbol"])

plotDF$Colour[(plotDF$log2FC < logFC_threshold & plotDF$log2FC > (-1*logFC_threshold)) | plotDF$adjPvalue > p_threshold] <- "unchanged"
plotDF$Colour[plotDF$log2FC < (-1*logFC_threshold) & plotDF$adjPvalue <= p_threshold] <- "downregulated"
plotDF$Colour[plotDF$log2FC >= logFC_threshold & plotDF$adjPvalue <= p_threshold] <- "upregulated"


selGenes <- c("PRSS35", "PLAT", "DACT1", "MEDAG", "CDKN2B", "TFPI",
              "SPOCD1", "CYP1B1", "OLFM2", "COL21A1")
# selGenes <- c(topTable_sel$hgnc_symbol[topTable_sel$logFC < 0][1:5],
#               topTable_sel$hgnc_symbol[topTable_sel$logFC > 0][1:5])

set.seed(123)
volcano <- ggplot2::ggplot() +
  ggplot2:: geom_point(data = plotDF[plotDF$Colour == "unchanged",],
                       ggplot2::aes(x = log2FC, y = logP, color = "Unchanged",
                                    alpha = logP, size = logP)) +
  ggplot2::geom_point(data = plotDF[plotDF$Colour == "downregulated",],
                      ggplot2::aes(x = log2FC, y = logP, color = "Downregulated",
                                   alpha = logP, size = logP)) +
  ggplot2::geom_point(data = plotDF[plotDF$Colour == "upregulated",],
                      ggplot2::aes(x = log2FC, y = logP, color = "Upregulated",
                                   alpha = logP, size = logP),) +
  #ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "lightgrey") +
  ggrepel::geom_text_repel(data = plotDF[plotDF$Name %in% selGenes,],
                           ggplot2::aes(x = log2FC, y = logP,
                                        label = Name),
                           min.segment.length = 0.1, force_pull = 0.2) +
  ggplot2::scale_color_manual(values = setNames(c(unchanged_color, down_color, up_color),
                                                c("Unchanged", "Downregulated", "Upregulated"))) +
  ggplot2::xlab(expression(log[2]~"FC")) +
  ggplot2::ylab(expression(-log[10]~"P value")) +
  ggplot2::xlim(-10,10) +
  ggplot2::ylim(0,23) +
  ggtitle("PAR1 KO: S vs C") +
  guides(size = "none", alpha = "none") +
  ggplot2::labs(color = NULL) + 
  ggplot2::scale_size_continuous(range = c(1,3)) +
  ggplot2::scale_alpha_continuous(range = c(0.5,1)) +
  ggplot2::theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

volcano

ggsave(volcano, width = 6, height = 5, file = "Output/Volcano_KO_alt.png")



#==============================================================================#
# Combine results into venn diagram
#==============================================================================#

# Collect up and downregulated genes
up_WT <- topTable_WT$ID[(topTable_WT$FDR < 0.05) & (topTable_WT$logFC > 0)]
down_WT <- topTable_WT$ID[(topTable_WT$FDR < 0.05) & (topTable_WT$logFC < 0)]
up_KO <- topTable_KO$ID[(topTable_KO$FDR < 0.05) & (topTable_KO$logFC > 0)]
down_KO <- topTable_KO$ID[(topTable_KO$FDR < 0.05) & (topTable_KO$logFC < 0)]

genesets <- list(`Upregulated\nin PAR1 WT` = up_WT, 
                 `Downregulated\nin PAR1 WT`=down_WT, 
                 `Upregulated\nin PAR1 KO`=up_KO,
                 `Downregulated\nin PAR1 KO`=down_KO)


# Make venn diagram
venn <- ggVennDiagram::ggVennDiagram(genesets, label_alpha = 0, label = "count") +
  coord_cartesian(xlim = c(0,1)) +
  theme(legend.position = "none") +
  ggplot2::scale_fill_gradient(low="white",high = "#CB181D")

ggsave(venn, width = 7, height = 5, file = "Output/Volcano_venn.png")

