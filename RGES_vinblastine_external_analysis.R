
#library(qvalue)
library("ROCR")

cancer <- "BRCA"
landmark <- 0
#CMAP score output
#GSE69845_drug_MCF7.RData
#GSE69850 Hepg2
output_path <- paste(cancer, "/lincs_score_", "GSE69850", ".csv", sep="")

##############
load("GSE69845_drug_MCF7.RData")
#load("GSE69850_drug_Hepg2.RData")

cmap_signatures <- drug_sig

lincs_drug_activity <- read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity <- unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",		"cell_line")))

gene.list <- cmap_signatures$GeneID

sig.ids <- names(cmap_signatures)[-c(1,2)]

#load  genes

load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC
res$GeneID <- sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol <- sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature <- subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf )
dz_signature <- subset(dz_signature, GeneID %in% gene.list ) #& GeneID %in% rownames(lincs_signatures))
dz_signature$up_down <- "up"
dz_signature$up_down[dz_signature$log2FoldChange<0] <- "down"


dz_genes_up <- subset(dz_signature,up_down=="up",select="GeneID")
dz_genes_down <- subset(dz_signature,up_down=="down",select="GeneID")

#only choose the top 100 genes
max_gene_size <- 150
if (nrow(dz_genes_up)> max_gene_size){
  dz_genes_up <- data.frame(GeneID= dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down)> max_gene_size){
  dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}


dz_cmap_scores <- NULL
dz_spearmans <- NULL
IQRs <- NULL
count <- 0
for (exp_id in sig.ids) {
  count <- count + 1
  print(paste("Computing score for disease against cmap_experiment_id =",count))
  sig <-  cmap_signatures[, exp_id]
  cmap_exp_signature <- data.frame(gene.list,  rank(-1 * sig, ties.method="random"))    
  colnames(cmap_exp_signature) <- c("ids","rank")
  dz_cmap_scores <- c(dz_cmap_scores, cmap_score_new(dz_genes_up,dz_genes_down,cmap_exp_signature))
  
 # compute spearman correlation and pearson correlation
  dz_sig_drug_sig <- merge(dz_signature, data.frame(GeneID=gene.list, expr= sig), by="GeneID")
  dz_spearmans <- c(dz_spearmans, cor(dz_sig_drug_sig$log2FoldChange, dz_sig_drug_sig$expr, method="spearman"))
  IQRs <- c(IQRs, IQR(sig))
}

results <- data.frame(id = sig.ids, RGES = dz_cmap_scores, spearman = dz_spearmans, iqr = IQRs)

results <- results[order(results$RGES),]
write.csv(results, output_path)

##########
#correlated to IC50
lincs_drug_activity <- aggregate(standard_value ~ pert_iname, lincs_drug_activity,median)
results$id <- tolower(results$id)
drug_activity_rges <- merge(results, lincs_drug_activity, by.x="id", by.y="pert_iname")

drug_sig <- cmap_signatures[,c("GeneID", "Vinblastine")]
