
require(xlsx)
require(pheatmap)
all_cell_lines_landmark = read.xlsx("raw/performance.xlsx", sheetName = "all_cell_lines_landmark")
all_cell_lines_landmark_t = data.frame(t(all_cell_lines_landmark[, -1]))
colnames(all_cell_lines_landmark_t) = all_cell_lines_landmark[,1]
all_cell_lines_landmark_t$cell_line = "all"
all_cell_lines_landmark_t$gene = "landmark"
all_cell_lines_landmark_t$measure = rownames(all_cell_lines_landmark_t)

cancer_specific_landmark = read.xlsx("raw/performance.xlsx", sheetName = "same_lineage_landmark")
cancer_specific_landmark_t = data.frame(t(cancer_specific_landmark[, -1]))
colnames(cancer_specific_landmark_t) = cancer_specific_landmark[,1]
cancer_specific_landmark_t$cell_line = "cancer"
cancer_specific_landmark_t$gene = "landmark"
cancer_specific_landmark_t$measure = rownames(cancer_specific_landmark_t)

all_cell_lines_all_genes = read.xlsx("raw/performance.xlsx", sheetName = "all_cell_lines_imputed")
all_cell_lines_all_genes_t = data.frame(t(all_cell_lines_all_genes[, -1]))
colnames(all_cell_lines_all_genes_t) = all_cell_lines_all_genes[,1]
all_cell_lines_all_genes_t$cell_line = "all"
all_cell_lines_all_genes_t$gene = "all"
all_cell_lines_all_genes_t$measure = rownames(all_cell_lines_all_genes_t)

cancer_specific_all_genes = read.xlsx("raw/performance.xlsx", sheetName = "same_lineage_imputed")
cancer_specific_all_genes_t = data.frame(t(cancer_specific_all_genes[, -1]))
colnames(cancer_specific_all_genes_t) = cancer_specific_all_genes[,1]
cancer_specific_all_genes_t$cell_line = "cancer"
cancer_specific_all_genes_t$gene = "all"
cancer_specific_all_genes_t$measure = rownames(cancer_specific_all_genes_t)

performance = rbind(all_cell_lines_landmark_t) #, cancer_specific_landmark_t,all_cell_lines_all_genes_t ,cancer_specific_all_genes_t)

performance$ave = apply(performance[, c(1,2,3)], 1, mean)
performance = performance[order(performance$ave), ]

performance_cmap = performance[c("Cmap.min", "Cmap.median", "Cmap.mean", "sRGES_other_lineage", "cell_sRGES_same_lineage", "sRGES_all"), ]
rownames(performance_cmap) = c("Best RGES", "Median RGES", "Mean RGES", "sRGES, other lineages", "sRGES, same lineage", "sRGES, all cell lines")

my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = 50)
pheatmap(performance_cmap[, c(1,2,3)], color = my_palette, cluster_rows=F, cluster_cols=F, 
         cellwidth=25,  cellheight=25, display_numbers=T, file="fig/performance.pdf")
