#correlate drug sensitivity with the potency to reserse disease gene expression (RGES)
#to reduce the computation time, only the drugs with sensitivty data were examined

#library(qvalue)
library("ROCR")
library("lsa")

cancer <- "BRCA"
landmark <- 1
#CMAP score output
#output_path <- paste(cancer, "/lincs_score_", landmark, ".csv", sep="")

if (cell_line_selected == "HT29"){
  cell_line_selected_chembl <- "HT-29"
}else if (cell_line_selected == "MCF7"){
  cell_line_selected_chembl <- "MCF7"
}else if (cell_line_selected == "HEPG2"){
  cell_line_selected_chembl <- "HepG2"
}

##############
load("raw/lincs/lincs_signatures_cmpd_landmark.RData")
lincs_sig_info <- read.csv("raw/lincs/lincs_sig_info.csv", stringsAsFactors=F)
lincs_drug_activity <- read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity <- unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",		"cell_line")))


if (landmark ==1){
  gene.list <- rownames(lincs_signatures)
}else{
  gene.list <- get.gene.list(con)  
}
sig.ids <- lincs_sig_info$id[lincs_sig_info$pert_type == "trt_cp" & lincs_sig_info$is_gold == 1 & lincs_sig_info$id >0 & tolower(lincs_sig_info$pert_iname) %in% tolower(lincs_drug_activity$pert_iname)]


#load  genes

load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC
res$GeneID <- sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol <- sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature <- subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf )
dz_signature <- subset(dz_signature, GeneID %in% gene.list )
dz_signature$up_down <- "up"
dz_signature$up_down[dz_signature$log2FoldChange<0] <- "down"

dz_genes_up <- subset(dz_signature,up_down == "up",select = "GeneID")
dz_genes_down <- subset(dz_signature,up_down == "down",select = "GeneID")

#only choose  top  genes
max_gene_size <- 150
if (nrow(dz_genes_up)> max_gene_size){
  dz_genes_up <- data.frame(GeneID= dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down)> max_gene_size){
  dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}


dz_cmap_scores <- NULL
dz_pearsons <- NULL
dz_spearmans <- NULL
IQRs <- NULL
cosines <- NULL
count <- 0
for (exp_id in sig.ids) {
  count <- count + 1
  print(paste("Computing score for disease against cmap_experiment_id =",count))
  if (landmark ==1){
    sig <-  lincs_signatures[, as.character(exp_id)]
    cmap_exp_signature <- data.frame(gene.list,  rank(-1 * sig, ties.method="random"))    
  }else{
    sig <-  get.instance.sig(exp_id, con)
    cmap_exp_signature <- cbind(gene.list,  rank(-1 * sig, ties.method="random"))    
  }
  colnames(cmap_exp_signature) <- c("ids","rank")
  dz_cmap_scores <- c(dz_cmap_scores, cmap_score_new(dz_genes_up,dz_genes_down,cmap_exp_signature))
  
}


results <- data.frame(id = sig.ids, RGES = dz_cmap_scores)

results <- merge(results, lincs_sig_info, by = "id")
results <- results[order(results$RGES),]



