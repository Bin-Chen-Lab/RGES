
#library(qvalue)
library("ROCR")

cancer = "BRCA"
landmark = 0
#CMAP score output
#GSE69845_drug_MCF7.RData
#GSE69850 Hepg2
output_path <- paste(cancer, "/lincs_score_", "GSE69850", ".csv", sep="")

##############
cmap_score_new <- function(sig_up, sig_down, drug_signature) {
  #the old function does not support the input list with either all up genes or all down genes, this new function attempts to addess this, not fully validated
  #Note. I think that creating the anonymous functions in each iteration of the sapply's below is slowing things down. Predefine them eventually.
  #we also tweak a bit: whenever the sign of ks_up/ks_down, we minus the two scores.
  num_genes = nrow(drug_signature)
  ks_up <- 0
  ks_down <- 0
  connectivity_score <- 0
  
  # I think we are re-ranking because the GeneID mapping changed the original rank range
  drug_signature[,"rank"] <- rank(drug_signature[,"rank"])
  
  # Merge the drug signature with the disease signature by GeneID. This becomes the V(j) from the algorithm description
  up_tags_rank <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
  down_tags_rank <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)
  
  up_tags_position <- sort(up_tags_rank$rank)
  down_tags_position <- sort(down_tags_rank$rank)
  
  num_tags_up <- length(up_tags_position)
  num_tags_down <- length(down_tags_position)
  
  # 
  if(num_tags_up > 1) {
    a_up <- 0
    b_up <- 0
    
    a_up <- max(sapply(1:num_tags_up,function(j) {
      j/num_tags_up - up_tags_position[j]/num_genes
    }))
    b_up <- max(sapply(1:num_tags_up,function(j) {
      up_tags_position[j]/num_genes - (j-1)/num_tags_up
    }))
    
    if(a_up > b_up) {
      ks_up <- a_up
    } else {
      ks_up <- -b_up
    }
  }else{
    ks_up <- 0
  }
  
  if (num_tags_down > 1){
    
    a_down <- 0
    b_down <- 0
    
    a_down <- max(sapply(1:num_tags_down,function(j) {
      j/num_tags_down - down_tags_position[j]/num_genes
    }))
    b_down <- max(sapply(1:num_tags_down,function(j) {
      down_tags_position[j]/num_genes - (j-1)/num_tags_down
    }))
    
    if(a_down > b_down) {
      ks_down <- a_down
    } else {
      ks_down <- -b_down
    }
  }else{
    ks_down <- 0
  }
  
  if (ks_up == 0 & ks_down != 0){ #only down gene inputed
    connectivity_score <- -ks_down 
  }else if (ks_up !=0 & ks_down == 0){ #only up gene inputed
    connectivity_score <- ks_up
  }else if (sum(sign(c(ks_down,ks_up))) == 0) {
    connectivity_score <- ks_up - ks_down # different signs
  }else{
    connectivity_score <- ks_up - ks_down
  }
  
  return(connectivity_score)
}

##############
load("GSE69845_drug_MCF7.RData")
#load("GSE69850_drug_Hepg2.RData")

cmap_signatures = drug_sig

lincs_drug_activity = read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity = unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",		"cell_line")))

gene.list = cmap_signatures$GeneID

sig.ids = names(cmap_signatures)[-c(1,2)]


#load  genes

load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC
res$GeneID = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf )
dz_signature = subset(dz_signature, GeneID %in% gene.list ) #& GeneID %in% rownames(lincs_signatures))
dz_signature$up_down = "up"
dz_signature$up_down[dz_signature$log2FoldChange<0] = "down"


dz_genes_up <- subset(dz_signature,up_down=="up",select="GeneID")
dz_genes_down <- subset(dz_signature,up_down=="down",select="GeneID")

#only choose the top 100 genes
max_gene_size = 150
if (nrow(dz_genes_up)> max_gene_size){
  dz_genes_up <- data.frame(GeneID= dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down)> max_gene_size){
  dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}


dz_cmap_scores = NULL
dz_spearmans = NULL
IQRs = NULL
count = 0
for (exp_id in sig.ids) {
  count = count + 1
  print(paste("Computing score for disease against cmap_experiment_id =",count))
  sig <-  cmap_signatures[, exp_id]
  cmap_exp_signature <- data.frame(gene.list,  rank(-1 * sig, ties.method="random"))    
  colnames(cmap_exp_signature) <- c("ids","rank")
  dz_cmap_scores = c(dz_cmap_scores, cmap_score_new(dz_genes_up,dz_genes_down,cmap_exp_signature))
  
 # compute spearman correlation and pearson correlation
  dz_sig_drug_sig = merge(dz_signature, data.frame(GeneID=gene.list, expr= sig), by="GeneID")
  dz_spearmans = c(dz_spearmans, cor(dz_sig_drug_sig$log2FoldChange, dz_sig_drug_sig$expr, method="spearman"))
  IQRs = c(IQRs, IQR(sig))
}

results = data.frame(id = sig.ids, RGES = dz_cmap_scores, spearman = dz_spearmans, iqr = IQRs)

results = results[order(results$RGES),]
write.csv(results, output_path)

##########
#correlated to IC50
lincs_drug_activity = aggregate(standard_value ~ pert_iname, lincs_drug_activity,median)
results$id = tolower(results$id)
drug_activity_rges = merge(results, lincs_drug_activity, by.x="id", by.y="pert_iname")

drug_sig = cmap_signatures[,c("GeneID", "Vinblastine")]
