# Cabo das: Limma Analysis of phosphoproteomic data 

# The purpose of this script is to analyze mass spec phosphoproteomic data,
# using the R package Limma to compare a full model (including terms for
# indvidual drugs as well as an interaction effect) to a reduced model
# that only accounts for individual drug effects.

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

filtro <- duplicated(pST$Phosphopeptide)
pST.use <- pST[!filtro, 10:17]
rownames(pST.use) <- pST$Phosphopeptide[!filtro]

# Fitting linear model to the pAll data
fit.pY   <- lmFit(pY.use, design)
fit.pY   <- eBayes(fit.pY)
fit.pST  <- lmFit(pST.use, design)
fit.pST  <- eBayes(fit.pST)

# Getting results, where logFC is for the cabo-das interaction effect
res.pY  <- topTable(fit.pY, coef = "Cabo:Das", sort.by = "none", resort.by = "M", n=Inf, p.value = 0.01)
res.pST <- topTable(fit.pST, coef = "Cabo:Das", sort.by = "none", resort.by = "M", n=Inf, p.value = 0.01)

res.pY %>% dim
res.pST %>% dim

res.pY  <- res.pY %>% 
  mutate(Phosphopeptide = rownames(res.pY)) %>% 
  dplyr::select(Phosphopeptide, logFC)

res.pST <- res.pST %>%
  mutate(Phosphopeptide = rownames(res.pST)) %>% 
  dplyr::select(Phosphopeptide, logFC)

