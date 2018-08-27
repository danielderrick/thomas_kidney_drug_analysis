# PROJECT : Cabo-Das
# TITLE   : Differential Expression analysis
#
# The purpose of this script is to import RNAseq counts, format them for
# DESeq2, identify differentially expressed genes, and export the processed
# data for downstream analysis.


# Setup
############################################################################
library(DESeq2)
library(tidyverse)
library(stringr)
library(tximport)

source("../../mddRNA/scripts/functions/DEfunctions.R") # handy


dir <- "~/data/thomas/cabo_das/ACHN/rnaseq/"
counts <- read.table(paste(dir, "counts/ACHN_cl_DAS.CABO.CYT.MK_rawcounts.txt",
                           sep = ""),
                     header = TRUE, 
                     row.names = 1, 
                     colClasses = c("character", rep("integer", 15)))


# Making DESeq Data Set
############################################################################

# Removing cyt_mk samples
counts <- counts[, -grep("cyt_mk", colnames(counts))]

# Making coldata
coldata <- colnames(counts) %>% str_split_fixed(., "[.]", 2)
colnames(coldata) <- c("drug", "id")
rownames(coldata) <- colnames(counts)

coldata <- 
  data.frame(coldata) %>% 
  mutate(sample = rownames(coldata)) %>% 
  mutate(cabo = case_when(
    grepl("cabo", drug) ~ "1",
    !grepl("cabo", drug) ~ "0")) %>% 
  mutate(das = case_when(
    grepl("das", drug) ~ "1",
    !grepl("das", drug) ~ "0")) %>% 
  dplyr::select(sample, das, cabo)

# Full model is effect of das + cabo + cabo:das
dds <- DESeqDataSetFromMatrix(counts, coldata, ~ das*cabo)

ACHN.dascabo.dds <- DESeq(dds, test = "LRT", reduced = ~ das + cabo)

out.dir <- "processed_data"

if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}


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

save(ACHN.dascabo.dds, 
     file = sprintf("%s/%s", out.dir, "ACHN_dascabo_dds.Rdata")

     