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

pAll <- rbind(pY.use, pST.use)

# Fitting linear model to the pAll data
fit.pAll.int  <- lmFit(pAll, design)
fit.pAll.int  <- eBayes(fit.pAll.int)

res.pAll.forheat.int  <- topTable(fit.pAll.int, 
                            coef = "Cabo:Das", sort.by = "none", 
                            resort.by = "M", n=Inf, p.value = 0.01,
                            lfc = 2.5)
pep.int <- rownames(res.pAll.forheat.int)
###############################################################################

meta <- data.frame(colnames(pY)[10:17])
colnames(meta) <- "sample"
meta$Cabo <- c(1, 1, 0, 0, 0, 0, 0, 0)
meta$Das  <- c(0, 0, 1, 1, 0, 0, 0, 0)
meta$CD   <- c(0, 0, 0, 0, 0, 0, 1, 1)
design <- model.matrix(~ Cabo + Das + CD, meta)

# Fitting linear model to the pAll data
fit.pAll.cd <- lmFit(pAll, design)
fit.pAll.cd <- eBayes(fit.pAll.cd)
res.pAll.forheat.cd  <- topTable(fit.pAll.cd, c("CD"), p.value = .05, lfc = 2.5, 
                                 number=Inf)

pep.cd <- rownames(res.pAll.forheat.cd)
# toPlot <- base::intersect(pep.int, pep.cd)
toPlot <- pep.cd


pY.w <- 
  pY %>% 
  mutate(type = "pY") %>% 
  select(-Sequence7, -Sequence10,
         -UniProt_ID, -Description, -Function.Phosphoresidue..phosphosite.org.,
         -Putative.Upstream.Kinases.Phosphatases.Binding.Domains)

pST.w <- 
  pST %>% 
  mutate(type = "pST") %>% 
  select(-Sequence7, -Sequence10,
       -UniProt_ID, -Description, -Function.Phosphoresidue..phosphosite.org.,
       -Putative.Upstream.Kinases.Phosphatases.Binding.Domains)

pData <- bind_rows(pY.w, pST.w) %>% 
  filter(Phosphopeptide %in% toPlot)

genefilt <- str_split(pData$Gene_Name, "; ") %>% lapply(., unique) %>% lapply(., length) %>% unlist
pData <- pData[genefilt == 1, ]

pData <- 
  pData %>%
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
  mutate(Phosphoresidue = sits) %>% 
  mutate(new_rowname = paste(Gene_Name, Phosphoresidue, sep = " "))

pData.toPlot <- 
  pData %>% 
  filter(!duplicated(new_rowname))

rownames(pData.toPlot) <- pData.toPlot$new_rowname
pData.toPlot <- 
  pData.toPlot %>% 
  dplyr::select(contains("Intensity"))

colnames(pData.toPlot) <- c("Veh A", "Veh B",
                            "Das A", "Das B",
                            "Cabo A", "Cabo B",
                            "Cabo + Das A",
                            "Cabo + Das B")


Heatmap(pData.toPlot, split = pData.toPlot$type,
        name = "intensity",
        cluster_columns = FALSE,
        col = colorRamp2(c(-2, 0, 2), c("red", "white", "blue")))
