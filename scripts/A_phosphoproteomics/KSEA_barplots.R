library(tidyverse)
library(cowplot)
library(forcats)

KSEA.pST <- read.table("/Users/derrickd/workspace/thomas/cabo_das/KSEA_new/KSEA_pST_og.txt",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
KSEA.pY  <- read.table("/Users/derrickd/workspace/thomas/cabo_das/KSEA_new/KSEA_pY_og.txt",
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

colnames(KSEA.pST)
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
  mutate(Motif = fct_inorder(Motif)) 

KSEA.names.tofix <- 
  KSEA.pST %>% 
  bind_rows(KSEA.pY) %>% 
  filter(FDR < 0.01) %>% 
  filter(No..hits > 5) %>% 
  filter(!duplicated(Motif)) %>% 
  select(Motif, Normalized_KS.Score, FDR, type) %>% 
  arrange(Normalized_KS.Score) %>% 
  mutate(Motif = fct_inorder(Motif)) %>% 
  pull(Motif) %>% 
  as.character()

new_names <- 
  data.frame(Motif = str_remove_all(as.character(KSEA.names.tofix), "\\(.+\\)"),
  HPRD = grepl("HPRD", as.character(KSEA.names.tofix)))

new_names <- 
  new_names %>% 
  mutate(HPRD = case_when(
    HPRD == TRUE ~ "(HPRD)",
    HPRD == FALSE ~ "")) %>% 
  mutate(Motif = paste(Motif, HPRD, sep = "")) %>% 
  pull(Motif)

KSEA.toPlot$Motif <- str_replace_all(KSEA.toPlot$Motif, "kinase", "kin.")
KSEA.toPlot$Motif <- str_replace_all(KSEA.toPlot$Motif, "substrate", "substr.")
KSEA.toPlot$Motif <- str_remove_all(KSEA.toPlot$Motif, "motif")
KSEA.toPlot$Motif <- str_remove_all(KSEA.toPlot$Motif, "Lumped")

KSEA.toPlot <-
  KSEA.toPlot %>% 
  arrange(Normalized_KS.Score) %>% 
  filter(!duplicated(Motif)) %>% 
  mutate(Motif = fct_inorder(Motif))
  
  
KSEA.toPlot <-
  KSEA.toPlot %>%
  filter(!duplicated(str_trim(as.character(KSEA.toPlot$Motif))))

pdf("~/Desktop/jul23/KSEA_figure.pdf", height = 10, width = 10)
ggplot(KSEA.toPlot, aes(x = Motif, y = Normalized_KS.Score, fill = type))  + geom_col(size = 4) + 
  theme(axis.text.y = element_text(size = 10)) +
  ggtitle("Kinase-Substrate Enrichment Analysis") + ylab("Normalized KS Score") + xlab("Kinase-Substrate Set") + coord_flip()
dev.off()
