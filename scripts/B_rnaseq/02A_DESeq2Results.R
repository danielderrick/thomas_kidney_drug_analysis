# PROJECT : Cabo-Das
# TITLE   : DESeq2 Results to Enrichr, RNK

###############################################################################
# Setup 
library(DESeq2)
library(tidyverse)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(eulerr)

load("processed_data/ACHN_dascabo_dds.Rdata")
source("../../mddRNA/scripts/functions/DEfunctions.R") # handy
###############################################################################

# Getting results
res <- results(ACHN.dascabo.dds)

# Adding symbols and filtering
res <- addSymbols2Res(data.frame(res))

res.filt <- res %>% 
  filter(abs(log2FoldChange) > .5) %>% 
  filter(padj < .05)

res.meta <-
  res.filt %>% 
  dplyr::select(ens, symbol)

achn.counts      <- assay(normTransform(ACHN.dascabo.dds))
achn.counts.norm <- t(scale(t(achn.counts), scale = FALSE))

save(achn.counts, achn.counts.norm,
     res, res.filt, res.meta, 
     ACHN.dascabo.dds,
     file = "processed_data/ACHN_rnaseq.Rdata")
