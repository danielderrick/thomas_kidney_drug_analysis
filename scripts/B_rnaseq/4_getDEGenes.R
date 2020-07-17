# PROJECT : Cabo-Das
# TITLE   : Differential Expression analysis
#
# The purpose of this script is to import RNAseq counts, format them for
# DESeq2 and identify DE genes for each treatment.
#
# Setup
############################################################################
library(eulerr)
library(DESeq2)
library(tidyverse)
library(stringr)
library(biomaRt)
library(GGally)


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

getG <- function(results, direction = NULL) {
  if (is.null(direction)) {
    print("Which direction? ('up' or 'down')")
  } else if (direction == "up") {
    x <-
      results %>% 
      data.frame %>% 
      addAnnotation(tx2g) %>% 
      filter(padj <= padj_threshold) %>% 
      filter(log2FoldChange >= log2FC_threshold) %>% 
      pull(ensembl_gene)
    x
  } else if (direction == "down") {
    x <-
      results %>% 
      data.frame %>% 
      addAnnotation(tx2g) %>% 
      filter(padj <= padj_threshold) %>% 
      filter(log2FoldChange <= -log2FC_threshold) %>% 
      pull(ensembl_gene)
    x
  } else {
    print("Which direction? ('up' or 'down')")
  }
}


# Preparing table for mapping gene identifiers
############################################################################
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = 'uswest.ensembl.org')

tx2g <- getBM(attributes = c("ensembl_gene_id",
                             "hgnc_symbol", 
                             "entrezgene_id"), 
              mart = mart) %>% 
  dplyr::rename(ensembl_gene = ensembl_gene_id,
                entrez_gene  = entrezgene_id)
tx2g[tx2g$hgnc_symbol == "", 2] <- NA

# Making DESeq Data Set
############################################################################

# Removing cyt-mk-treated samples - not examined in this analysis
counts <- counts[, -grep("cyt_mk", colnames(counts))]

# Making coldata
meta <- 
  colnames(counts) %>% 
  str_split_fixed(., "[.]", 2)

dimnames(meta) <- list(rows = colnames(counts), 
                       columns = c("drug", "id"))

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
  dplyr::select(sample, drug, das, cabo)

############################################################################
# Creating DESeq Data Set - first, with single factor for drug
dds <- DESeqDataSetFromMatrix(counts, coldata, ~ drug)

# Running DESeq comparing full model (das + cabo + das*cabo) vs reduced model (das + cabo)
dds <- DESeq(dds)

dc <- lfcShrink(dds, contrast = c("drug", "das_cabo", "veh"))
d  <- lfcShrink(dds, contrast = c("drug", "das", "veh"))
c  <- lfcShrink(dds, contrast = c("drug", "cabo", "veh"))

padj_threshold <- 0.01
log2FC_threshold <- 0.5

up <- list(
  "Cabozantinib" = getG(c, "up"),
  "Dasatinib"  = getG(d, "up"),
  "Combination" = getG(dc, "up")
  )

down <- list(
  "Cabozantinib" = getG(c, "down"),
  "Dasatinib"  = getG(d, "down"),
  "Combination" = getG(dc, "down")
)

# Venn diagrams
pdf(sprintf("../venn_diagram_folder/venn_nocolor_genes_logfc%s.pdf", log2FC_threshold), height = 4.5, width = 4.5)
euler(up) %>% 
  plot(quantities = TRUE,
       fills = c("white", "white", "white"),
       main = sprintf("P < %s; log2FC > %s", padj_threshold, log2FC_threshold))

euler(down) %>% 
  plot(quantities = TRUE,
       fills = c("white", "white", "white"),
       main = sprintf("P < %s; log2FC < -%s", padj_threshold, log2FC_threshold))
dev.off()

###############################################################################
dcTJ <-
  dc %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  mutate(DE = padj <= padj_threshold &
           abs(log2FoldChange) >= log2FC_threshold) %>% 
  dplyr::rename(L2FC = log2FoldChange) %>% 
  dplyr::select(ensembl_gene_id, contains("L2FC"), contains("DE")) %>% 
  dplyr::rename(DE_Combo = DE,
                L2FC_Combo = L2FC)

dTJ <-
  d %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  mutate(DE = padj <= padj_threshold &
           abs(log2FoldChange) >= log2FC_threshold) %>% 
  dplyr::rename(L2FC = log2FoldChange) %>% 
  dplyr::select(ensembl_gene_id, contains("L2FC"), contains("DE")) %>% 
  mutate(treatment = "Das")

cTJ <-
  c %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  mutate(DE = padj <= padj_threshold &
           abs(log2FoldChange) >= log2FC_threshold) %>% 
  dplyr::rename(L2FC = log2FoldChange) %>% 
  dplyr::select(ensembl_gene_id, contains("L2FC"), contains("DE")) %>% 
  mutate(treatment = "Cabo")

cTJ %>% 
  bind_rows(dTJ) %>% 
  pivot_wider(names_from = treatment, values_from = c(DE, L2FC)) %>% 
  mutate(L2FC_Added = L2FC_Das + L2FC_Cabo) %>% 
  pivot_longer(-ensembl_gene_id, 
               names_to = c(".value", "treatment"),
               names_sep = "_") %>% 
  left_join(dcTJ) %>%
  ggplot(aes(L2FC_Combo, L2FC, color = DE_Combo)) +
  geom_point(shape = 16, alpha = 0.25) +
  geom_smooth(color = "black", method = "lm") +
  facet_wrap(~treatment) +
  stat_regline_equation()

l2fcMat <- 
  cTJ %>% 
  bind_rows(dTJ) %>% 
  bind_rows(dplyr::rename(dcTJ, L2FC = L2FC_Combo,
                          DE = DE_Combo) %>% 
              mutate(treatment = "Combination")) %>% 
  dplyr::select(-DE) %>% 
  pivot_wider(names_from = treatment, values_from = L2FC)  %>% 
  filter(!is.na(Cabo),
         !is.na(Das),
         !is.na(Combination)) %>% 
  column_to_rownames("ensembl_gene_id")

cor(l2fcMat)^2

pdf("figuresRemade/Treatment_L2FC_Scatterplots.pdf", height = 6, width = 6)
ggpairs(l2fcMat,
        # mapping = aes(alpha = .1),
        lower=list(continuous= wrap("smooth", alpha = .1, size = .5)),
        diag=list(continuous="bar"), 
        upper=list(continuous = wrap("cor", method = "pearson"))
        )
dev.off()


###############################################################################
# Numbers for text:

# Cabo genes:
list("up" = up$Cabozantinib, 
     "down" = down$Cabozantinib, 
     "both" = c(up$Cabozantinib, down$Cabozantinib)) %>% 
  lapply(length)

# Das genes:
list("up" = up$Dasatinib, 
     "down" = down$Dasatinib, 
     "both" = c(up$Dasatinib, down$Dasatinib)) %>% 
  lapply(length)

# Combination genes
list("up" = up$Combination, 
     "down" = down$Combination, 
     "both" = c(up$Combination, down$Combination)) %>% 
  lapply(length)

###############################################################################
# Recreating DESeq DataSet with two-factor design
# Full model is effect of das + cabo + cabo:das
dds    <- DESeqDataSetFromMatrix(counts, coldata, ~ das*cabo)

# Running DESeq comparing full model (das + cabo + das*cabo) vs reduced model (das + cabo)
dds <- DESeq(dds, 
             test = "LRT", 
             reduced = ~ das + cabo)

sigInteraction <-
  results(dds) %>% 
  data.frame %>% 
  rownames_to_column("ensembl_gene") %>% 
  filter(padj <= 0.01)

sigInteraction %>% dim
