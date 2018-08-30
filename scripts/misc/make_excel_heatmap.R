library(gdata)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

mats <- list("ACHN\n% Viability"          = read.xls("misc/drug_matrix.xlsx", sheet = 1),
             "ACHN\n% Growth Inhibition"  = read.xls("misc/drug_matrix.xlsx", sheet = 2),
             "CAKI\n% Viability"          = read.xls("misc/drug_matrix.xlsx", sheet = 3),
             "CAKI\n% Growth Inhibition"  = read.xls("misc/drug_matrix.xlsx", sheet = 4),
             "786\n% Viability"          = read.xls("misc/drug_matrix.xlsx", sheet = 5),
             "786\n% Growth Inhibition"  = read.xls("misc/drug_matrix.xlsx", sheet = 6),
             "SN12C\n% Viability"         = read.xls("misc/drug_matrix.xlsx", sheet = 7),
             "SN12C\n% Growth Inhibition" = read.xls("misc/drug_matrix.xlsx", sheet = 8))

mats_fixed <- 
  lapply(mats, function(x) {
  colnames(x) <- as.character(x[1, ])
  rownames(x) <- x[, 1]
  x <- x[-1, -1]
  x
})

cols <- c("indianred", "white", "dodgerblue1")
scales.matchExcel <- list("achn_v"  = colorRamp2(c(.15, .67, 1.14), cols),
                          "achn_gi" = colorRamp2(c(0,   .32,  .85), rev(cols)),
                          "caki_v"  = colorRamp2(c(.02, .50, 1.01), cols),
                          "caki_gi"  = colorRamp2(c(.00, .61, .98), rev(cols)),
                          "786_v"  = colorRamp2(c(.02, .19, 1.01), cols),
                          "786_gi"  = colorRamp2(c(0.0, .78, .98), rev(cols)),
                          "sn12c_v"  = colorRamp2(c(10.6, 38.46, 99.7), cols),
                          "sn12c_gi"  = colorRamp2(c(0, .62, .89), rev(cols)))



i <- 1L
pdf("excel_heatmaps.pdf", width = 8, height = 8)
par(mfrow = c(2,2))
mapply(function(x, i) {
  Heatmap(x,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               name = names(mats_fixed)[i],
               col = scales.matchExcel[[i]],
               row_title = "[Dasatanib] (µM) ",
               column_title_side = "top",
               column_names_side = "top",
               column_title = "[Cabozantinib] (µM)",
               row_names_side = "left",
               row_title_side = "left")
}, x = mats_fixed[c(2,4,6,8)], i = c(2, 4, 6, 8))

dev.off()
