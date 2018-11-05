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
meta$Cabo <- c(0, 0, 0, 0, 1, 1, 0, 0)
meta$Das  <- c(0, 0, 1, 1, 0, 0, 0, 0)
meta$CD   <- c(0, 0, 0, 0, 0, 0, 1, 1)

design <- model.matrix(~ Cabo + Das + CD, meta)
pY.use <- pY[, 10:17]
pY.use$type <- "pY"
rownames(pY.use) <- pY$Phosphopeptide

filtro <- duplicated(pST$Phosphopeptide)
pST.use <- pST[!filtro, 10:17]
pST.use$type <- "pST"
rownames(pST.use) <- pST$Phosphopeptide[!filtro]

pAll <- rbind(pY.use, pST.use)
key <- setNames(pAll$type, rownames(pAll))
pAll <- pAll[, -9]

# Fitting linear model to the pAll data
fit.pAll.int  <- lmFit(pAll, design)
fit.pAll.int  <- eBayes(fit.pAll.int)

res.pAll.C  <- topTable(fit.pAll.int,
                        coef = "Cabo", sort.by = "none",
                        resort.by = "M", n=Inf, p.value = 0.01,
                        lfc = 2)
rall.c <- rownames(res.pAll.C)
split.c <- key[rall.c]

res.pAll.D  <- topTable(fit.pAll.int,
                        coef = "Das", sort.by = "none",
                        resort.by = "M", n=Inf, p.value = 0.01,
                        lfc = 2)
rall.d <- rownames(res.pAll.D)
split.d <- key[rall.d]

res.pAll.CD  <- topTable(fit.pAll.int,
                        coef = "CD", sort.by = "none",
                        resort.by = "M", n=Inf, p.value = 0.01,
                        lfc = 2)
rall.cd <- rownames(res.pAll.CD)
split.cd <- key[rall.cd]

colnames(pAll) <- c("Veh A", "Veh B",
                    "Das A", "Das B",
                    "Cabo A", "Cabo B",
                    "Cabo + Das A",
                    "Cabo + Das B")

dir <- c("figures/heatmaps/single_agent/phospho/")

if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

pdf("figures/heatmaps/single_agent/phospho/das_phospho_heatmap.pdf", height = 8, width = 6)
Heatmap(pAll[rall.d, ],
        split = split.d,
        show_row_names = FALSE,
        col = colorRamp2(c(-2.25, 0, 2.25), c("blue", "white", "red")),
        name = "Intensity",
        column_title = "Dasatinib\nlog2FC > 2; q < 0.01",
        cluster_columns = FALSE)
dev.off()

pdf("figures/heatmaps/single_agent/phospho/cabo_phospho_heatmap.pdf", height = 8, width = 6)
Heatmap(pAll[rall.c, ],
        split = split.c,
        show_row_names = FALSE,
        col = colorRamp2(c(-2.25, 0, 2.25), c("blue", "white", "red")),
        name = "Intensity",
        column_title = "Cabozantinib\nlog2FC > 2; q < 0.01",
        cluster_columns = FALSE)
dev.off()

pdf("figures/heatmaps/single_agent/phospho/cabodas_phospho_heatmap.pdf", height = 8, width = 6)
Heatmap(pAll[rall.cd, ],
        split = split.cd,
        show_row_names = FALSE,
        col = colorRamp2(c(-2.25, 0, 2.25), c("blue", "white", "red")),
        name = "Intensity",
        column_title = "Cabozantinib+Dasatinib\nlog2FC > 2; q < 0.01",
        cluster_columns = FALSE)
dev.off()
