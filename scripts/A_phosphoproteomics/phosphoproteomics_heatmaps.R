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
  mutate(Diff = CD - Ctrl) %>% 
  mutate(type = "pY") %>% 
  select(-Phosphopeptide, -Sequence7, -Sequence10,
         -UniProt_ID, -Description, -Function.Phosphoresidue..phosphosite.org.,
         -Putative.Upstream.Kinases.Phosphatases.Binding.Domains)

pST.w <- 
  pST %>% 
  mutate(CD = (Intensity.CD.4A + Intensity.CD.4B)/2) %>%
  mutate(Ctrl = (Intensity.CTL.1A + Intensity.CTL.1B)/2) %>% 
  mutate(Diff = CD - Ctrl) %>% 
  mutate(type = "pST") %>% 
  select(-Phosphopeptide, -Sequence7, -Sequence10,
       -UniProt_ID, -Description, -Function.Phosphoresidue..phosphosite.org.,
       -Putative.Upstream.Kinases.Phosphatases.Binding.Domains)

pData <- bind_rows(pY.w, pST.w)

genefilt <- str_split(pData$Gene_Name, "; ") %>% lapply(., unique) %>% lapply(., length) %>% unlist
pData <- pData[genefilt == 1, ]

pData <- 
  pData %>% 
  filter(abs(Diff) > 1.5) %>% 
  mutate(Phosphoresidue = str_remove_all(Phosphoresidue, "[;,]"),
         Gene_Name = str_extract(Gene_Name, "[:alnum:]+"))

sits <- str_split(pData$Phosphoresidue, " ")
sits <- sapply(sits, unique)

sits <- sapply(sits, function(x) {
  if (!is.na(x[4])) {
    x[3] <- paste(x[3], "*", sep = "")
    x <- paste(x[1:3], sep = ",", collapse = ", ")
  } else {
    x <- paste(x, sep = ",", collapse = ", ")
  }
  x
})
pData <- 
  pData %>% 

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
