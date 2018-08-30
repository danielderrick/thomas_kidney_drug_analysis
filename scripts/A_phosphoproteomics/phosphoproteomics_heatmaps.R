# Cabo das: Phosphoproteomic data heatmaps

# The purpose of this script is to create heat maps for the phosphoproteomic 
# data for cabo, das, and combination-treated ACHN cells.

###############################################################################
# Setup
library(gdata)
library(stringr)
library(reshape2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

pY <- read.xls("~/data/thomas/cabo_das/ACHN/phosphoproteomics/intensity/201712_GT_pY_invitro_ExcelOutput.xlsx",
                sheet = 4)

pST <- read.xls("~/data/thomas/cabo_das/ACHN/phosphoproteomics/intensity/201712_GT_pST_invitro_ExcelOutput.xlsx",
                sheet = 4)
###############################################################################

# Getting difference between combination and control measurements
pY.w <- 
  pY %>% 
  mutate(CD = (Intensity.CD.4A + Intensity.CD.4B)/2) %>%
  mutate(Ctrl = (Intensity.CTL.1A + Intensity.CTL.1B)/2) %>% 
  mutate(Diff = CD - Ctrl)

pST.w <- 
  pST %>% 
  mutate(CD = (Intensity.CD.4A + Intensity.CD.4B)/2) %>%
  mutate(Ctrl = (Intensity.CTL.1A + Intensity.CTL.1B)/2) %>% 
  mutate(Diff = CD - Ctrl)

pY.filt <-
  pY.w %>% 
  filter(abs(Diff) > 1.5) %>%
  mutate(type = "Y") %>% 
  dplyr:: select(Phosphopeptide, type,
         Intensity.CTL.1A, Intensity.CTL.1B, 
         Intensity.Cabo.2A, Intensity.Cabo.2B,
         Intensity.Das.3A, Intensity.Das.3B,
         Intensity.CD.4A, Intensity.CD.4B)
rownames(pY.filt) <- pY.filt$Phosphopeptide
pY.plot <-
  pY.filt %>% 
  dplyr::select(-type, -Phosphopeptide)

pST.filt <-
  pST.w %>% 
  filter(abs(Diff) > 1.5) %>% 
  mutate(type = "ST") %>% 
  filter(!duplicated(Phosphopeptide)) %>% 
  dplyr::select(Phosphopeptide, type,
                 Intensity.CTL.1A, Intensity.CTL.1B, 
                 Intensity.Cabo.2A, Intensity.Cabo.2B,
                 Intensity.Das.3A, Intensity.Das.3B,
                 Intensity.CD.4A, Intensity.CD.4B)

pData <- bind_rows(pY.filt, pST.filt) %>% 
  filter(!duplicated(Phosphopeptide))
pData$type <- factor(pData$type,
                     levels = c("ST", "Y"),
                     labels = c("pST", "pY"))

pData.toPlot <- pData
rownames(pData.toPlot) <- pData$Phosphopeptide
pData.toPlot <- 
  pData.toPlot %>% 
  dplyr::select(-Phosphopeptide, -type)

colnames(pData.toPlot) <- c("Veh A", "Veh B", "Cabo A", "Cabo B", "Das A", "Das B", "Cabo + Das A", "Cabo + Das B")

pdf("plots/heatmaps/phosphodata/phosphodata_heatmap_combined.pdf")
Heatmap(pData.toPlot, split = pData$type,
        show_row_names = FALSE,
        name = "intensity",
        col = colorRamp2(c(-2, 0, 2), c("red", "white", "blue")))
dev.off()
