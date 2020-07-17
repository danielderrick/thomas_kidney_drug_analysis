# Cabo Das: Generate upstream (kinase) TieDie input
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
KSEA.out <- read.csv("/Users/derrickd/data/thomas/cabo_das/ACHN/phosphoproteomics/KSEA/interaction/interaction_pAll.csv")

outDir <- "processed_data/TieDie/TieDie_input"
###############################################################################

# Making sure that pAll is all character/numeric vectors
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

temp.pY  <- topTable(fit.pY, coef = "Cabo:Das", p.value = .05, n =Inf)
temp.pST <- topTable(fit.pST, coef = "Cabo:Das", p.value = .01, n =Inf)

temp.pY.join <-
  temp.pY %>% 
  mutate(Phosphopeptide = rownames(temp.pY)) %>% 
  dplyr::select(Phosphopeptide, logFC) %>% 
  inner_join(., pY, by = "Phosphopeptide")

temp.pST.join <-
  temp.pST %>% 
  mutate(Phosphopeptide = rownames(temp.pST)) %>% 
  dplyr::select(Phosphopeptide, logFC) %>% 
  inner_join(., pST, by = "Phosphopeptide")

temp.pY.join$Function.Phosphoresidue..phosphosite.org. %>% 
  levels %>% 
  grep("phosphor", ., value = TRUE)

filtro.pY  <- grep("phosphorylat", 
                   as.character(temp.pY.join$Function.Phosphoresidue..phosphosite.org.))
filtro.pST <- grep("phosphorylat", 
                   as.character(temp.pST.join$Function.Phosphoresidue..phosphosite.org.))

tpAll.j <- rbind({temp.pY.join[filtro.pY, ] %>% dplyr::select(Gene_Name, logFC)},
                 {temp.pST.join[filtro.pST, ] %>% dplyr::select(Gene_Name, logFC)})

tpAll.j$Gene_Name <- str_extract(tpAll.j$Gene_Name, "[:alnum:]+")
tpAll.j           <- tpAll.j[-c(16,17), ]

direct <- tpAll.j

####################################### Adding KSEA output

KSEA.out <- KSEA.out[!is.na(KSEA.out$gene_name), ]
KSEA.out <- KSEA.out[!duplicated(KSEA.out$gene_name), ]
colnames(KSEA.out) <- c("gene", "logfc")
colnames(direct)   <- c("gene", "logfc")

all.upstream <- bind_rows(KSEA.out, direct)

filtro <- duplicated(all.upstream$gene)
all.upstream <- all.upstream[!filtro, ]

upstream.df <- data.frame(cbind(all.upstream$gene,
                 sign(all.upstream$logfc),
                 abs(all.upstream$logfc)),
                stringsAsFactors = FALSE)

upstream.df[upstream.df$X2 == -1, 2] <- "-"
upstream.df[upstream.df$X2 == 1, 2] <- "+"
upstream.df[, 3] <- round(as.numeric(upstream.df[, 3]), 4)
upstream.df <- upstream.df %>% dplyr::select(X1, X3, X2)

####################################### Writing to File

if(!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

write.table(upstream.df,
            file = sprintf("%s/upstream.input", outDir),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
