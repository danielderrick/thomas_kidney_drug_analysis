# PROJECT : Cabo-Das
# TITLE   : Generating Master Regulator NES Scores

##############    WHAT    #################################################
# The purpose of this script is to calculate master regulator normalized
# enrichment scores for ACHN cells treated with dasatanib and cabozantinib.

##############    WHY    ##################################################
# We are interested in master regulators that change in activity after
# combination drug treatment. MRs implicated will be used as TieDie input.

##############    HOW    ##################################################
# Counts are imported. VIPER and a papillary RCC specific interactome are
# used to generate normalized enrichment scores.

##############    INPUT    ################################################
# ACHN ± Das/Cabo DESeq data               ACHN_4Drug_DESeq2.Rdata
# Papillary RCC interactome                regulonkirp

library(viper)
library(aracne.networks)
library(tidyverse)
library(org.Hs.eg.db)

data("regulonkirp")
load("processed_data/ACHN_rnaseq.Rdata")

###########################################################################

vp.sig <- 
  res %>% 
  filter(!is.na(padj)) %>% 
  dplyr::select(ens, log2FoldChange) %>% 
  mutate(entrez = mapIds(org.Hs.eg.db, ens, "ENTREZID", "ENSEMBL")) %>% 
  filter(!is.na(entrez)) %>% 
  filter(!is.na(log2FoldChange)) %>% 
  dplyr::select(entrez, log2FoldChange)

vp.sig <- setNames(vp.sig$log2FoldChange, nm = vp.sig$entrez)

vp.res <- msviper(vp.sig, regulonkirp)

tmp <- unique(c(names(vp.res$regulon), 
                rownames(vp.res$signature)))
annot <- mapIds(org.Hs.eg.db, tmp, "SYMBOL", "ENTREZID")
vp.res <- msviperAnnot(vp.res, annot)

vp.res.syn <- msviperCombinatorial(vp.res, regulators = 400)
vp.res.syn <- msviperSynergy(vp.res.syn)
plot(vp.res.syn, mrs = 50)
vp.res.syn$es$synergy

pdf("interaction_mr_analysis.pdf", height = 7, width = 10)
plot(vp.res, mrs = 30)
dev.off()


filtro <- (abs(vp.res$es$nes) > 5)

downstream.toWrite <- vp.res$es$nes[filtro, drop = FALSE]

downstream.toWrite <- data.frame(
  cbind(names(downstream.toWrite),
      sign(downstream.toWrite),
      abs(downstream.toWrite)),
  stringsAsFactors = FALSE)

downstream.toWrite[downstream.toWrite$X2 == -1, 2] <- "-"
downstream.toWrite[downstream.toWrite$X2 == 1, 2] <- "+"
downstream.toWrite <- dplyr::select(downstream.toWrite, X1, X3, X2)
downstream.toWrite$X3 <- round(as.numeric(downstream.toWrite$X3), digits = 4)

write.table(downstream.toWrite,
            file = "~/Desktop/downstream.input",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
