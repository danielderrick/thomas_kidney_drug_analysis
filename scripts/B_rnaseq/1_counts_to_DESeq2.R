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
library(biomaRt)

dir <- "~/data/thomas/cabo_das/ACHN/rnaseq/"
counts <- read.table(paste(dir, "counts/ACHN_cl_DAS.CABO.CYT.MK_rawcounts.txt",
                           sep = ""),
                     header = TRUE, 
                     row.names = 1, 
                     colClasses = c("character", rep("integer", 15)))

symbolizeResults <- function(x, key) {
  # a function to take a results table with ensembl ids as rownames and
  # convert to a data frame with a "hgnc_symbol" column
  x <-
    data.frame(x)
  x <- 
    x %>%
    mutate(ensembl_gene = rownames(x)) %>% 
    inner_join(key, by = "ensembl_gene") %>% 
    filter(!duplicated(hgnc_symbol)) %>% 
    filter(!is.na(hgnc_symbol)) %>% 
    filter(!is.na(padj)) %>% 
    filter(!is.na(log2FoldChange)) %>% 
    dplyr::select(hgnc_symbol, everything()) %>% 
    dplyr::select(-ensembl_gene, entrez_gene)
  x
}

entrezResults <- function(x, key) {
  # a function to take a results table with ensembl ids as rownames and 
  # convert to a data frame with a "hgnc_symbol" column
  x <-
    data.frame(x)
  x <- 
    x %>% 
    mutate("ensembl_gene" = rownames(x)) %>% 
    inner_join(key, by = "ensembl_gene") %>%
    filter(!duplicated(entrez_gene)) %>%
    filter(!is.na(entrez_gene)) %>%
    filter(!is.na(log2FoldChange)) %>%
    dplyr::select(entrez_gene, everything()) %>%
    dplyr::select(-ensembl_gene, hgnc_symbol)
  x
}

addAnnotation <- function(x, tx2g) {
  x <- data.frame(x) %>% 
    mutate(ensembl_gene = rownames(x))
  tx2g <- 
    tx2g %>% 
    dplyr::select(ensembl_gene, entrez_gene, hgnc_symbol) %>% 
    filter(!duplicated(ensembl_gene))
  x <- 
    left_join(x, tx2g, by = "ensembl_gene") %>% 
    dplyr::select(ensembl_gene, entrez_gene, hgnc_symbol, everything())
  return(x)
}

writeCounts <- function(x, filename) {
  write.table(x,
              file = filename,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE) 
}

# Preparing table for mapping gene identifiers
############################################################################
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = 'www.ensembl.org')

tx2g <- getBM(attributes = c("ensembl_gene_id",
                             "hgnc_symbol", 
                             "entrezgene"), 
              mart = mart) %>% 
  dplyr::rename(ensembl_gene = ensembl_gene_id,
                entrez_gene  = entrezgene)
tx2g[tx2g$hgnc_symbol == "", 2] <- NA

# Making DESeq Data Set
############################################################################

# Removing cyt-mk-treated samples - not examined in this analysis
counts <- counts[, -grep("cyt_mk", colnames(counts))]

# Making coldata
meta <- 
  colnames(counts) %>% 
  str_split_fixed(., "[.]", 2) %>% 
  data.frame()

dimnames(meta) <- list(rows = colnames(counts), 
                       columns = c("drug", "id"))
meta$drug <- fct_shift(meta$drug, 3)
# rearranging coldata into a two factor design - samples are either 
# yes or no (0 or 1) for das and cabo
coldata <- 
  data.frame(meta) %>% 
  mutate(sample = rownames(meta)) %>% 
  mutate(cabo = case_when(
    grepl("cabo", drug) ~ "1",
    !grepl("cabo", drug) ~ "0")) %>% 
  mutate(das = case_when(
    grepl("das", drug) ~ "1",
    !grepl("das", drug) ~ "0")) %>% 
  mutate(das = as.factor(das),
         cabo = as.factor(cabo)) %>% 
  dplyr::select(sample, das, cabo)

# Full model is effect of das + cabo + cabo:das
dds    <- DESeqDataSetFromMatrix(counts, coldata, ~ das*cabo)
dds2   <- DESeqDataSetFromMatrix(counts, meta, ~ drug)

# Running DESeq comparing full model (das + cabo + das*cabo) vs reduced model (das + cabo)
dds <- DESeq(dds, 
             test = "LRT", 
             reduced = ~ das + cabo)

# Getting counts ##############################################################
achn.rnaseq <- vector("list", length = 2)

# counts normalized for lib size differences
achn.rnaseq[[1]] <- assay(normTransform(dds)) 
# counts normalized for lib size differences and mean-centered
achn.rnaseq[[2]] <- t(scale(t(achn.rnaseq[[1]]), scale = FALSE))

names(achn.rnaseq) <- c("norm", "mean_centered")

achn.rnaseq <- lapply(achn.rnaseq, function(x) {
  x <- addAnnotation(x, tx2g = tx2g)
  x
})

# Getting results #############################################################

# getting results - log2FoldChange is for interaction effect of das + cabo 
res <- vector("list", length = 3L)

res[[1]] <- results(dds)
res[[2]] <- symbolizeResults(res[[1]], key = tx2g)
res[[3]] <- entrezResults(res[[1]], key = tx2g)

names(res) <- c("ensembl", "gene_symbol", "entrez")

res.filt <- lapply(res, function(x) {
  x <- 
    data.frame(x) %>% 
    filter(padj < .01)
  x
})

out.dir <- "processed_data"

if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}

save(dds,
     dds2,
     res,
     res.filt,
     achn.rnaseq,
     file = sprintf("%s/%s",
                    out.dir, "rnaseq_data.Rdata"))

writeCounts(achn.rnaseq$norm,
            sprintf("%s/%s",
                    out.dir, "norm_counts.tsv"))

writeCounts(achn.rnaseq$norm,
            sprintf("%s/%s",
                    out.dir, "norm_counts_meancent.tsv"))
