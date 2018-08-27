counts.tp1 <- counts[helper.tab$ens, ]
rownames(counts.tp1) <- helper.tab$symbol

# Heatmap for the ~ genes
pdf("interaction_heatmap.pdf", width = 8, height = 12)
Heatmap(counts.tp1,
        name = "mean-centered\nlog2(counts+1)",
        row_names_gp = gpar(fontsize = 8))
dev.off()

# writing log2 fold changes for GSEA
write.table(interaction,
            file = "interaction.rnk",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)

# Writing threshold-passing genes for Enrichr
res.filt %>% 
  filter(!is.na(symbol)) %>% 
  filter(!is.na(log2FoldChange)) %>% 
  filter(padj < .05) %>% 
  filter(log2FoldChange > 0) %>% 
  pull(symbol) %>% 
  writeLines(., "forEnrichrUp.txt")


res.filt %>% 
  filter(!is.na(symbol)) %>% 
  filter(!is.na(log2FoldChange)) %>% 
  filter(padj < .05) %>% 
  filter(log2FoldChange < 0) %>% 
  pull(symbol) %>% 
  writeLines(., "forEnrichrDown.txt")

res.filt %>% pull(log2FoldChange) %>% hist


res.filt %>% arrange(desc(abs(log2FoldChange))) %>% filter(!is.na(symbol)) %>%  filter(!is.na(padj)) %>% head(., n = 10)
