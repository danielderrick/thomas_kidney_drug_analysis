library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

load("processed_data/rnaseq_data.Rdata")
rm(dds2, res, res.filt)

filtro <- duplicated(achn.rnaseq$mean_centered$hgnc_symbol) | is.na(achn.rnaseq$mean_centered$hgnc_symbol)
achn.rnaseq$mean_centered <- achn.rnaseq$mean_centered[!filtro, ]
rownames(achn.rnaseq$mean_centered) <- achn.rnaseq$mean_centered$hgnc_symbol
achn.rnaseq$mean_centered <- achn.rnaseq$mean_centered[, -c(1:3)]

# Functions for annotating results
############################################################################
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

makePalette <- function(data, squish_range = c(0.001, .999)) {
  # The purpose of this function is to make a ColorRamp2 palette by
  # "squishing" extreme values, then using the squished values and the palette
  # from Mark to produce a palette for ComplexHeatmap
  pal_choice <- "RdBu"
  
  squish <- function(mat, lower_thresh = .005, upper_thresh = .995) {
    # The purpose of this function is to "squish" a matrix.
    # Values in the matrix that are more extreme than input quantiles
    # are "squished" down to the quantile value.
    range <- quantile(unlist(mat), c(lower_thresh, upper_thresh))
    
    if (range[1] > -1) {
      range[1] <- -1
    }
    
    if (range[2] < 1) {
      range[2] <- 1
    }
    
    mat <- apply(mat, c(1,2), function(x) {
      if (x > range[2]) {
        x <- range[2]
      } else if (x < range[1]) {
        x <- range[1]
      }
      return(x)
    })
    return(mat)
  }
  
  pal <- brewer.pal(n=11,pal_choice)
  rc1 <- colorRampPalette(colors = c(pal[1], pal[2]), space="Lab")(11)
  for(i in 2:10){
    tmp <- colorRampPalette(colors = c(pal[i], pal[i+1]), space = "Lab")(10)
    rc1 <- c(rc1,tmp)
  }
  
  if (!is_null(squish_range)) {
    data <- squish(data, 
                   lower_thresh = squish_range[1], 
                   upper_thresh = squish_range[2])
  }
  
  rb1 <- seq(min(data), 0, length.out=50+1)
  rb2 <- seq(0, max(data), length.out=50+1)[-1]
  rampbreaks <- c(rb1, rb2)
  palette <- colorRamp2(breaks = rampbreaks,
                        colors = rev(rc1))
}

# Preparing table for mapping gene identifiers
############################################################################
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = 'www.ensembl.org')

tx2g <- getBM(attributes = c("ensembl_gene_id",
                             "hgnc_symbol"), 
              mart = mart) %>% 
  dplyr::rename(ensembl_gene = ensembl_gene_id)
tx2g[tx2g$hgnc_symbol == "", 2] <- NA

# Preparing table for mapping gene identifiers
############################################################################
dds <- DESeq(dds, test = "LRT", reduced = ~ das + cabo)

res <- results(dds)

res.f <-
  res %>% 
  data.frame() %>% 
  mutate(ensembl_gene = rownames(res)) %>% 
  filter(padj < 0.05) %>% 
  filter(abs(log2FoldChange) > 0.5) %>% 
  left_join(tx2g, by = "ensembl_gene")

res.f.symbols <- res.f %>% pull(hgnc_symbol)
res.f.symbols <- res.f.symbols[!is.na(res.f.symbols)]

res.f.symbols.up <- res.f %>%
  filter(log2FoldChange > 0) %>% 
  pull(hgnc_symbol)

res.f.symbols.down <- res.f %>%
  filter(log2FoldChange < 0) %>% 
  pull(hgnc_symbol)

writeLines(res.f.symbols.up, "~/Desktop/dascabo_int_genes_up.txt")
writeLines(res.f.symbols.down, "~/Desktop/dascabo_int_genes_down.txt")

dir <- "figures/heatmaps/interaction/rnaseq/"
if(!dir.exists(dir)) {
  dir.create(dir)
}

pdf("figures/heatmaps/interaction/rnaseq/interaction_rnaseq_heatmap.pdf", height = 11, width = 6)
Heatmap(achn.rnaseq$mean_centered[res.f.symbols, ],
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 7),
        name = "mean-centered\nlog2(counts + 1)",
        col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
        column_title = "Full vs reduced model\nlog2FC > 0.5")
dev.off()

length(res.f.symbols)

