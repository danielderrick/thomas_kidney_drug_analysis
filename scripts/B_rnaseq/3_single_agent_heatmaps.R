library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

load("processed_data/rnaseq_data.Rdata")
rm(dds, res, res.filt)
dds <- dds2

rownames(achn.rnaseq$mean_centered) <- achn.rnaseq$mean_centered$ensembl_gene
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
                             "hgnc_symbol", 
                             "entrezgene"), 
              mart = mart) %>% 
  dplyr::rename(ensembl_gene = ensembl_gene_id,
                entrez_gene  = entrezgene)
tx2g[tx2g$hgnc_symbol == "", 2] <- NA

# Preparing table for mapping gene identifiers
############################################################################
dds <- DESeq(dds)

thresh <- 1.5
resDC <- results(dds, contrast = c("drug", "veh", "das_cabo"),
                 lfcThreshold = thresh)
resD  <- results(dds, contrast = c("drug", "veh", "das"), 
                 lfcThreshold = thresh)
resC  <- results(dds, contrast = c("drug", "veh", "cabo"), 
                 lfcThreshold = thresh)

resList <- list("DC" = resDC,
                "D"  = resD,
                "C"  = resC)

genesList <- lapply(resList, function(x) {
  X <- data.frame(x) %>% 
    mutate(ens = rownames(x)) %>% 
    filter(padj < 0.05) %>% 
    pull(ens)
  X
})

dir <- c("figures/heatmaps/single_agent/lfc15")

if (!dir.exists(dir)) {
  dir.create(dir, recursive = TRUE)
}

pdf("figures/heatmaps/single_agent/lfc15/cabo_heatmap.pdf", height = 6, width = 6)
Heatmap(achn.rnaseq$mean_centered[genesList$C, ],
        show_row_names = FALSE,
        name = "mean-centered\nlog2(counts + 1)",
        col = colorRamp2(c(-2.25, 0, 2.25), c("blue", "white", "red")),
        column_title = "Cabo vs Vehicle\nlog2(FC) > 1.5")
dev.off()

pdf("figures/heatmaps/single_agent/lfc15/das_cabo_heatmap.pdf", height = 15, width = 6)
Heatmap(achn.rnaseq$mean_centered[genesList$DC, ],
        show_row_names = FALSE,
        name = "mean-centered\nlog2(counts + 1)",
        col = colorRamp2(c(-2.25, 0, 2.25), c("blue", "white", "red")),
        column_title = "Das+Cabo vs Vehicle\nlog2(FC) > 1.5")
dev.off()


cabo <- addAnnotation(achn.rnaseq$mean_centered[genesList$C, ],
                      tx2g)
filtro <- duplicated(cabo$hgnc_symbol) | is.na(cabo$hgnc_symbol)
cabo <- cabo[!filtro, ]
rownames(cabo) <- cabo$hgnc_symbol
cabo <- cabo[, -c(1:3)]

pdf("figures/heatmaps/single_agent/lfc15/cabo_heatmap_named.pdf", height = 6, width = 6)
Heatmap(cabo,
        row_names_gp = gpar(fontsize = 8),
        name = "mean-centered\nlog2(counts + 1)",
        col = colorRamp2(c(-2.25, 0, 2.25), c("blue", "white", "red")),
        column_title = "Cabo vs Vehicle\nlog2(FC) > 1.5")
dev.off()

dc <- addAnnotation(achn.rnaseq$mean_centered[genesList$DC, ],
                      tx2g)
filtro <- duplicated(dc$hgnc_symbol) | is.na(dc$hgnc_symbol)
dc <- dc[!filtro, ]
rownames(dc) <- dc$hgnc_symbol
dc <- dc[, -c(1:3)]

pdf("figures/heatmaps/single_agent/lfc15/das_cabo_heatmap_named.pdf", height = 15, width = 6)
Heatmap(dc,
        row_names_gp = gpar(fontsize = 5),
        name = "mean-centered\nlog2(counts + 1)",
        col = colorRamp2(c(-2.25, 0, 2.25), c("blue", "white", "red")),
        column_title = "Das+Cabo vs Vehicle\nlog2(FC) > 1.5")
dev.off()

