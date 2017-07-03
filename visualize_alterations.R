
library(ComplexHeatmap)
cancer <- "LIHC"
mat <- read.table(paste( "reversal_genes/", cancer, "_mutation.txt", sep=""), 
                 header = TRUE,stringsAsFactors=FALSE, sep = "\t")
mat[is.na(mat)] <- ""
rownames(mat) <- mat[, 1]
mat <- mat[, -1]
#mat=  mat[, -ncol(mat)]
mat <- t(as.matrix(mat))
rownames(mat) <- sapply(rownames(mat), function(x) unlist(strsplit(x, "\\."))[1])
mat[1:3, 1:3]

## -------------------------------------------------------------------------------------------------
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)

## -------------------------------------------------------------------------------------------------
col <- c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")

#column_order = sample_order

## ---- fig.width = 12, fig.height = 8--------------------------------------------------------------
pdf(paste("fig/oncoprint_", cancer, ".pdf", sep=""))

oncoPrint(mat, get_type = function(x) unlist(strsplit(strsplit(x, ";")[[1]], ":"))[1],
          alter_fun = alter_fun, col = col, row_order = NULL, remove_empty_columns = T, 
          column_title = "", show_row_barplot = T, 
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "HOMDEL", "MUT"), 
                                      labels = c("Amplification", "Deep deletion", "Mutation")))
dev.off()