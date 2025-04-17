# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(edgeR)
library(tidyverse)

# Set working directory
setwd("PATH/TO/DIRECTORY")

# Read and format count data
count_data <- read.delim("Data/CountMatrix.txt")
count_matrix = as.matrix(count_data[,-1])
rownames(count_matrix) = count_data[,1]

# Make meta data
metaData <- data.frame(SampleID = colnames(count_matrix),
                       Phenotype = factor(substring(unlist(lapply(str_split(colnames(count_matrix),"_"), function(x) x[[2]])),1,1)),
                       Genotype = factor(unlist(lapply(str_split(colnames(count_matrix),"_"), function(x) x[[1]]))),
                       Replicate = factor(substring(unlist(lapply(str_split(colnames(count_matrix),"_"), function(x) x[[2]])),2,2))
)

# Check if samples are in correct order
all(metaData$SampleID == colnames(count_matrix))

#*****************************************************************************#
# Normalize data (edgeR method)
#*****************************************************************************#

# Make DGE list object
pheno_geno <- factor(paste0(metaData$Phenotype, "_", metaData$Genotype))
y <- DGEList(counts = count_matrix,
             group = pheno_geno)

# Remove lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y)

# Transform the data
gxMatrix_norm <- log2(cpm(y) + 1)

# Save the data
save(gxMatrix_norm, file = "gxMatrix_norm.RData")


#*****************************************************************************#
# Perform statistical analysis
#*****************************************************************************#

# Make design matrix
design <- model.matrix(~0 + pheno_geno)
colnames(design) <- levels(pheno_geno)

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(y, design)

# Make Contrasts
all_contrasts <- makeContrasts(
  complex = (S_KO - C_KO) - (S_WT - C_WT), # Complex comparison
  SvsC_KO = S_KO - C_KO,                   # S vs C in KO
  SvsC_WT = S_WT - C_WT,                   # S vs C in WT
  levels = design
)

# Create top tables:

# 1. Complex comparison
test <- glmQLFTest(fit, contrast = all_contrasts[,"complex"])
topTable_complex <- topTags(test, n = Inf)$table
topTable_complex$ID <- rownames(topTable_complex)

# 2. S vs C in KO
test <- glmQLFTest(fit, contrast = all_contrasts[,"SvsC_KO"])
topTable_KO <- topTags(test, n = Inf)$table
topTable_KO$ID <- rownames(topTable_KO)

# 3. S vs C in WT
test <- glmQLFTest(fit, contrast = all_contrasts[,"SvsC_WT"])
topTable_WT <- topTags(test, n = Inf)$table
topTable_WT$ID <- rownames(topTable_WT)

# Save top table
save(topTable_complex, topTable_KO, topTable_WT, file = "topTable.RData")


#*****************************************************************************#
# Prepare data for PathVisio
#*****************************************************************************#

# Combine WT and KO top tables
outputTable <- inner_join(topTable_WT[,c(1,4,5,6)], topTable_KO[,c(1,4,5,6)],
                          by = c("ID" = "ID"))

# Set column names
colnames(outputTable) <- c("logFC_WT", "Pvalue_WT", "FDR_WT", "ID",
                           "logFC_KO", "Pvalue_KO", "FDR_KO")

# Re-order columns
outputTable <- outputTable[,c("ID", 
                              "logFC_WT", "Pvalue_WT", "FDR_WT",
                              "logFC_KO", "Pvalue_KO", "FDR_KO")]

# Add complex comparison to data frame
outputTable <- inner_join(outputTable, topTable_complex[,c(4,5,6)],
                          by = c("ID" = "ID"))

# Set column names
colnames(outputTable) <- c("ID","logFC_WT", "Pvalue_WT", "FDR_WT",
                           "logFC_KO", "Pvalue_KO", "FDR_KO",
                           "Pvalue_Complex", "FDR_Complex")

# Write table that can be used as input for PathVisio
write.table(outputTable, "PathVisioTable.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)
