
load("raw/lincs/lincs_signatures_cmpd_landmark.RData")

sig_info = read.csv("raw/lincs/lincs_sig_info.csv")
sig_info = sig_info[sig_info$pert_iname == "pyrvinium-pamoate" & sig_info$cell_id %in% c("HEPG2", "HUH7"),]

lincs_signatures_subset = lincs_signatures[, as.character(sig_info$id)]

load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #
res$GeneID = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf )

genes_affected = merge(dz_signature, data.frame(GeneID = names(lincs_signatures_subset), z_score = as.numeric(lincs_signatures_subset)), by = "GeneID")
genes_affected = genes_affected[order(genes_affected$z_score), ]

genes_selected =   genes_affected$symbol[abs(genes_affected$z_score) > 1.5 ]

wnt_genes = c("LRP6", "BIRC5", "CCND1", "AXIN1", "BCL2", "BCL2L1", "CTNNB1", "CSNK1G1", "CSNK1G2", "CSNK1G3", "CDKN1A")
load("~/Documents/stanford/tumor_cell_line/data/GSE62944/TumorFPKM.RData")
load("~/Documents/stanford/tumor_cell_line/data/GSE62944/clinvar.RData")
clinvar_subset = clinvar[clinvar$CancerType == "LIHC",]

plot(genes_affected$GeneID)
tumor_expr_case_select = t(TumorFPKM[rownames(TumorFPKM) %in% c( genes_selected), ])[rownames(clinvar_subset), ]
tumor_expr_case_wnt = t(TumorFPKM[rownames(TumorFPKM) %in% c( wnt_genes), ])[rownames(clinvar_subset), ]

cors  = cor(tumor_expr_case_select, tumor_expr_case_wnt)
cors
sort(apply(cors, 1, mean))

plot(log(tumor_expr_case_select[, "EGR1"]), log( tumor_expr_case_wnt[, "BIRC5"]))
pheatmap(cors, file = "fig/pyrvinium_co_expression.pdf")

sort(cors["BIRC5",])

cor.test(a[,1], a[,2])

genes_affected = merge(res, data.frame(GeneID = names(lincs_signatures_subset), z_score = as.numeric(lincs_signatures_subset)), by = "GeneID")
genes_affected = genes_affected[order(genes_affected$z_score), ]
genes_affected$seq = seq(1:nrow(genes_affected))
genes_affected$up_down = NA
genes_affected$up_down[genes_affected$log2FoldChange> 1.5 & genes_affected$padj < 0.001] =  "up"
genes_affected$up_down[genes_affected$log2FoldChange< -1.5 & genes_affected$padj < 0.001] =  "down"
genes_affected$point_size = ifelse(is.na(genes_affected$up_down), 1, 1.5 )
genes_affected$wnt = ifelse(genes_affected$symbol %in% wnt_genes, T, F)
#aes( color = up_down, size = point_size)
pdf("fig/pyrvinium_expr.pdf")
  ggplot(genes_affected, aes(x=seq, y=z_score)) +  theme_bw()  + 
  geom_point() + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  xlab("") + guides(shape=FALSE, size=FALSE) +
  ylab("Z score")
dev.off()

library("GSVA")
library("GSA")
ssgsea_matrix = data.frame(value = lincs_signatures_subset)
rownames(ssgsea_matrix) = names(lincs_signatures_subset)
ssgsea_matrix = as.matrix(ssgsea_matrix)
msigdb = GSA.read.gmt("~/Documents/stanford//breast/data/msigdb.v5.0.entrez.gmt")
geneSets = msigdb$genesets
names(geneSets) = msigdb$geneset.names
gsea_results = gsva(ssgsea_matrix, geneSets, method = "ssgsea")
gsea_results = gsea_results[order(gsea_results[,1]), ]

gsea_results[grep("WNT", names(gsea_results))][1:20]
gsea_results[grep("KEGG", names(gsea_results))][1:20]
gsea_results[grep("REACTOME", names(gsea_results))][1:20]
