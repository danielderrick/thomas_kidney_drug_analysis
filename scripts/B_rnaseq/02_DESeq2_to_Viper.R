# PROJECT : Cabo-Das
# TITLE   : Generating Master Regulator NES Scores
#
# The purpose of this script is to calculate master regulator normalized
# enrichment scores for ACHN cells treated with dasatanib and cabozantinib.


# Setup
###########################################################################
library(viper)
library(aracne.networks)
library(tidyverse)
library(biomaRt)

data("regulonkirp")
load("processed_data/rnaseq_data.Rdata")
###############################################################################

# Preparing table for mapping gene identifiers
###############################################################################
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

# Extracting interaction signature from DESeq2 results, running viper
###############################################################################
# getting signature as named vector of log2foldchanges (names are entrez ids)
vp.sig.prep <-
  res$entrez %>% 
  filter(!is.na(padj)) %>% 
  filter(!is.na(entrez_gene)) %>% 
  filter(!is.na(log2FoldChange)) %>% 
  dplyr::select(entrez_gene, log2FoldChange)
vp.sig <- setNames(vp.sig.prep$log2FoldChange, vp.sig.prep$entrez_gene)

# running viper using kirp regulon
vp.res <- msviper(vp.sig, regulonkirp)

# creating table to convert entrez ids to hgnc symbols, annotating results
annot.prep <- 
  tx2g %>% 
  filter(entrez_gene %in% unique(c(names(vp.res$regulon), 
                                   rownames(vp.res$signature)))) %>% 
  filter(!is.na(hgnc_symbol)) %>% 
  dplyr::select(entrez_gene, hgnc_symbol)
annot <- setNames(annot.prep$hgnc_symbol, 
                  as.character(annot.prep$entrez_gene))
vp.res.a   <- msviperAnnot(vp.res, annot)

out.dir <- "figures"

if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}

pdf(sprintf("%s/%s", out.dir, "interaction_mr_analysis.pdf", 
            height = 7, width = 10))
plot(vp.res.a, mrs = 30)
dev.off()


# Extracting master regulators with |NES| > 5 as input to TieDie
###############################################################################

filtro <- (abs(vp.res.a$es$nes) > 5)

downstream.toWrite <- vp.res$es$nes[filtro, drop = FALSE]

downstream.toWrite <- data.frame(
  cbind(names(downstream.toWrite),
      sign(downstream.toWrite),
      abs(downstream.toWrite)),
  stringsAsFactors = FALSE)

filtro <- downstream.toWrite$X2 == -1
downstream.toWrite[filtro,  2] <- "-"
downstream.toWrite[!filtro, 2]  <- "+"

downstream.toWrite <- 
  downstream.toWrite %>% 
  dplyr::select(X1, X3, X2) %>% 
  mutate(X3 = round(as.numeric(downstream.toWrite$X3), digits = 4))

tiedie.dir <- "tiedie/input"
if (!dir.exists(tiedie.dir)) {
  dir.create(tiedie.dir, recursive = TRUE)
}

write.table(downstream.toWrite,
            file = sprintf("%s/%s",
                           tiedie.dir, "downstream.input"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
