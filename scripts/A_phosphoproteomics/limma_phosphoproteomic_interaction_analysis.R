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
fit.pY  <- lmFit(pY.use, design)
fit.pY   <- eBayes(fit.pY)
fit.pST <- lmFit(pST.use, design)
fit.pST  <- eBayes(fit.pST)

# Getting results, where logFC is for the cabo-das interaction effect
res.pY  <- topTable(fit.pY, coef = "Cabo:Das", sort.by = "none", resort.by = "M", n=Inf)
res.pST <- topTable(fit.pST, coef = "Cabo:Das", sort.by = "none", resort.by = "M", n=Inf)

res.pY  <- res.pY %>% 
  mutate(Phosphopeptide = rownames(res.pY)) %>% 
  dplyr::select(Phosphopeptide, logFC)
s
res.pST <- res.pST %>%
  mutate(Phosphopeptide = rownames(res.pST)) %>% 
  dplyr::select(Phosphopeptide, logFC)

pY.toWrite <- 
  full_join(pY, res.pY, by = "Phosphopeptide") %>% 
  select(Phosphopeptide, logFC, everything())

pY.toWrite.sm <-
  pY.toWrite %>% 
  select(Phosphopeptide, logFC) %>% 
  arrange(desc(logFC))


pST.toWrite <- 
  full_join(pST, res.pST, by = "Phosphopeptide") %>% 
  select(Phosphopeptide, logFC, everything())

pST.toWrite.sm <-
  pST.toWrite %>% 
  select(Phosphopeptide, logFC) %>% 
  arrange(desc(logFC))

write.csv(pST.toWrite.sm, file = "pST_ranked_phosphopeptides.csv", quote = FALSE)
write.csv(pY.toWrite.sm, file = "pY_ranked_phosphopeptides.csv", quote = FALSE)


top <- res %>%
  mutate(phosphopeptide = rownames(res)) %>% 
  arrange(desc((logFC))) %>%
  pull(phosphopeptide) %>% 
  head(., n = 150)


toptab <- pAll.use[top, ]
filtro <- duplicated(pAll.map[top, 3]) | is.na(pAll.map[top, 3])

toptab <- toptab[!filtro, ]
rownames(toptab) <- pAll.map[rownames(toptab), 3]

pdf("toptab_desc.pdf", height = 15, width = 6)
Heatmap(toptab,
        row_names_gp = gpar(fontsize = 8),
        cluster_columns = FALSE)
dev.off()

top <- res %>%
  mutate(phosphopeptide = rownames(res)) %>% 
  arrange(((logFC))) %>%
  pull(phosphopeptide) %>% 
  head(., n = 150)


toptab <- pAll.use[top, ]
filtro <- duplicated(pAll.map[top, 3]) | is.na(pAll.map[top, 3])

toptab <- toptab[!filtro, ]
rownames(toptab) <- pAll.map[rownames(toptab), 3]

pdf("toptab_ascen.pdf", height = 15, width = 6)
Heatmap(toptab,
        row_names_gp = gpar(fontsize = 8),
        cluster_columns = FALSE)
dev.off()

top <- res %>%
  mutate(phosphopeptide = rownames(res)) %>% 
  arrange(desc(abs(logFC))) %>%
  pull(phosphopeptide) %>% 
  head(., n = 150)


toptab <- pAll.use[top, ]
filtro <- duplicated(pAll.map[top, 3]) | is.na(pAll.map[top, 3])

toptab <- toptab[!filtro, ]
rownames(toptab) <- pAll.map[rownames(toptab), 3]

pdf("toptab_desc_abs.pdf", height = 15, width = 6)
Heatmap(toptab,
        row_names_gp = gpar(fontsize = 8),
        cluster_columns = FALSE)
dev.off()


top <- res %>%
  mutate(phosphopeptide = rownames(res)) %>% 
  arrange((abs(logFC))) %>%
  pull(phosphopeptide) %>% 
  head(., n = 400)
top <- pep2

toptab <- pAll.use[top, ]
filtro <- duplicated(pAll.map[top, 3]) | is.na(pAll.map[top, 3])

toptab <- toptab[!filtro, ]
rownames(toptab) <- pAll.map[rownames(toptab), 3]

pdf("toptab_final.pdf", height = 8, width = 8)
Heatmap(toptab,
        row_names_gp = gpar(fontsize = 10),
        cluster_columns = FALSE,
        name = "intensity")
dev.off()

pep <- readLines("interesting_peptides2.txt")
pep2 <- sapply(pep, function(x) {
  str_split(x, pattern = " ")
})

pep2 <- unlist(pep2)
names(pep2) <- NULL

meta
colnames(pAll)
