library(tidyverse)
library(cowplot)
library(forcats)

theme_set(theme_cowplot())

KSEA.pST <- read.table("processed_data/KSEA_output/KSEA_pST_og.txt",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
KSEA.pY  <- read.table("processed_data/KSEA_output/KSEA_pY_og.txt",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# KSEA.pST <-
KSEA.pST <-
  KSEA.pST %>% 
  mutate(type = "pST")

KSEA.pY <-
  KSEA.pY %>% 
  mutate(type = "pY")

KSEA.toPlot <- 
  KSEA.pST %>% 
  bind_rows(KSEA.pY) %>% 
  filter(FDR < 0.01) %>% 
  filter(No..hits > 5) %>% 
  filter(!duplicated(Motif)) %>% 
  select(Motif, Normalized_KS.Score, FDR, type) %>% 
  arrange(Normalized_KS.Score) %>% 
  mutate(Motif = fct_inorder(Motif),
         Motif = str_remove_all(as.character(Motif), "\\(.+\\)")) %>% 
  mutate(HPRD = grepl("HPRD", as.character(KSEA.names.tofix))) %>% 
  mutate(HPRD = case_when(HPRD == TRUE ~ "(HPRD)", HPRD == FALSE ~ "")) %>% 
  mutate(Motif = paste(Motif, HPRD, sep = "")) %>% 
  mutate(Motif = str_replace_all(Motif, "kinase", "kin.")) %>% 
  mutate(Motif = str_replace_all(Motif, "substrate", "substr.")) %>% 
  mutate(Motif = str_remove_all(Motif, "motif")) %>% 
  mutate(Motif = str_remove_all(Motif, "Lumped")) %>% 
  mutate(Motif = str_remove_all(Motif, "Lumped")) %>% 
  mutate(Motif = str_trim(Motif)) %>% 
  arrange(Normalized_KS.Score) %>% 
  filter(!duplicated(Motif)) %>% 
  mutate(Motif = fct_inorder(Motif))

outDir <- "figuresRemade"

if(!dir.exists(outDir)) {dir.create(outDir)}

pdf(sprintf("%s/KSEA_interaction_barplot.pdf", outDir),
    height = 10, width = 10)
ggplot(KSEA.toPlot, aes(x = Motif, y = Normalized_KS.Score, fill = type))  + geom_col(size = 4) + 
  theme(axis.text.y = element_text(size = 10)) +
  ggtitle("Kinase-Substrate Enrichment Analysis") + ylab("Normalized KS Score") + xlab("Kinase-Substrate Set") + coord_flip()
dev.off()
