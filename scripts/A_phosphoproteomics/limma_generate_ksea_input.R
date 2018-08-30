# Cabo Das: Calculate Interaction log2FoldChanges for phosphoproteomic data
#
# The purpose of this script is to calculate log2FoldChanges for the
# phosphoproteomic data for the interaction effect. The output is then
# written to CSV files to be sent to Larry Cheng @ Rutgers for KSEA
# analysis.
#
###############################################################################
# Setup
library(gdata)
library(stringr)
library(reshape2)
library(tidyverse)
library(ComplexHeatmap)
library(limma)

pY <- read.xls("~/data/thomas/cabo_das/ACHN/phosphoproteomics/intensity/201712_GT_pY_invitro_ExcelOutput.xlsx",
               sheet = 1)

pST <- read.xls("~/data/thomas/cabo_das/ACHN/phosphoproteomics/intensity/201712_GT_pST_invitro_ExcelOutput.xlsx",
                sheet = 1)
###############################################################################

# Making sure that pAll is "all" character/numeric vectors
pY$Phosphopeptide <- as.character(pY$Phosphopeptide)
pY$Gene_Name <- as.character(pY$Gene_Name)

pST$Phosphopeptide <- as.character(pST$Phosphopeptide)
pST$Gene_Name <- as.character(pST$Gene_Name)

# Reordering columns for convenience
pY <- pY[, c(1:9, 14:17, 10:11, 12:13)]
pST <- pST[, c(1:9, 14:17, 10:11, 12:13)]

# Making metadata table for building model matrix
meta <- data.frame(colnames(pY)[10:17])
colnames(meta) <- "sample"
meta$Cabo <- c(0, 0, 0, 0, 1, 1, 1, 1)
meta$Das <- c(0, 0, 1, 1, 0, 0, 1, 1)
design <- model.matrix(~Cabo*Das, meta)

pY.use <- pY[, 10:17]
rownames(pY.use) <- pY$Phosphopeptide

pST.use <- pST[, 10:17]
rownames(pST.use) <- pST$Phosphopeptide

# Fitting linear model to the pAll data
fit.pY  <- lmFit(pY.use, design)
fit.pY   <- eBayes(fit.pY)
fit.pST <- lmFit(pST.use, design)
fit.pST  <- eBayes(fit.pST)

res.pY  <- topTable(fit.pY, coef = "Cabo:Das", sort.by = "none", resort.by = "M", n=Inf)
res.pST <- topTable(fit.pST, coef = "Cabo:Das", sort.by = "none", resort.by = "M", n=Inf)

res.pY.f  <- res.pY %>% 
  mutate(Phosphopeptide = rownames(res.pY)) %>% 
  dplyr::select(Phosphopeptide, logFC)

res.pST.f <- res.pST %>% 
  mutate(Phosphopeptide = rownames(res.pST)) %>% 
  dplyr::select(Phosphopeptide, logFC)

pY.toWrite <- 
  full_join(pY, res.pY.f, by = "Phosphopeptide") %>% 
  select(Phosphopeptide, logFC, everything())

pY.toWrite.sm <-
  pY.toWrite %>% 
  select(Phosphopeptide, logFC) %>% 
  arrange(desc(logFC))


pST.toWrite <- 
  full_join(pST, res.pST.f, by = "Phosphopeptide") %>% 
  select(Phosphopeptide, logFC, everything())

pST.toWrite.sm <-
  pST.toWrite %>% 
  select(Phosphopeptide, logFC) %>% 
  arrange(desc(logFC))

write.csv(pST.toWrite.sm, file = "processed_data/pST_ranked_phosphopeptides.csv", quote = FALSE)
write.csv(pY.toWrite.sm, file = "processed_data/pY_ranked_phosphopeptides.csv", quote = FALSE)
