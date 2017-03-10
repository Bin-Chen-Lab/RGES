#visualize HCC drugs
#
library(pheatmap)
library("gplots")
library(ConsensusClusterPlus)
library(metap)
library("RColorBrewer")
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)) #RdYlBu

selected_drug_names = c( "CGK-733",  "strophanthidin", "FCCP", "pyrvinium-pamoate")

cancer = "LIHC"
landmark = 1
fdr_cutoff = 0.25
#CMAP score output
output_path <- paste(cancer, "/all_lincs_score.csv", sep="")

lincs_drug_prediction = read.csv(output_path)
lincs_drug_prediction$RGES = lincs_drug_prediction$cmap_score

########
load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #
res$GeneID = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf )
dz_sig = subset(dz_signature, select=c("GeneID", "log2FoldChange"))

##########
#visualized gene reversed
#only pick the signatures from close to the median
cell_lines = read.csv(paste( "raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))
lincs_drug_prediction_subset = subset(lincs_drug_prediction, pert_iname %in% selected_drug_names, select=c("id", "RGES", "pert_iname", "pert_dose", "pert_time"))

###selecting median still sounds weird... let's keep all signatures 
drug_cmap_score = aggregate(RGES ~ pert_iname, lincs_drug_prediction_subset, median)
drug_instances_median  = merge(lincs_drug_prediction_subset, drug_cmap_score, by = c("pert_iname"))
drug_instances_median$diff = abs(drug_instances_median$RGES.x - drug_instances_median$RGES.y) #cmap_score.y is the median
drug_instances_min_diff = aggregate(diff ~ pert_iname, drug_instances_median, min)
drug_instances_select = merge(drug_instances_median, drug_instances_min_diff, by=c("pert_iname", "diff"))
drug_instances_select = drug_instances_select[!duplicated(drug_instances_select$pert_iname), ]
sig_id_selects = drug_instances_select$id

#sig_id_selects = lincs_drug_prediction_subset$id

load("raw/lincs/lincs_signatures_cmpd_landmark.RData")
#load(paste(cancer, "/", cancer, "_drug_signatures1.RData", sep=""))
#lincs_signatures = drug_signatures

drug_dz_signature = merge(dz_sig, data.frame(GeneID = rownames(lincs_signatures), lincs_signatures[, as.character(sig_id_selects)]),  by="GeneID", suffixes='')


#########################
###
#visualize the reversed gene expression
#reorder drug dz signatures
gene_ids = drug_dz_signature$GeneID
drug_dz_signature_rank = drug_dz_signature[,-1]
for (i in 1:ncol(drug_dz_signature_rank)){
  drug_dz_signature_rank[,i] = rank(-1 * drug_dz_signature_rank[,i] ) #highly expressed genes ranked on the top
}
gene_ids_rank <- gene_ids[order(drug_dz_signature_rank[,1])]
drug_dz_signature_rank <- drug_dz_signature_rank[order(drug_dz_signature_rank[,1]),] #order by disease expression

lincs_drug_activity = read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity = unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",  	"cell_line")))
lincs_drug_activity = aggregate(standard_value ~ pert_iname, lincs_drug_activity, median)
active_drug_names = lincs_drug_activity$pert_iname[lincs_drug_activity$standard_value <= 10000 & lincs_drug_activity$pert_iname %in% lincs_drug_prediction_subset$pert_iname]
inactive_drug_names = lincs_drug_activity$pert_iname[lincs_drug_activity$standard_value > 10000 & lincs_drug_activity$pert_iname %in% lincs_drug_prediction_subset$pert_iname]
#median_drug_names = lincs_drug_activity$pert_iname[lincs_drug_activity$standard_value > 10000 & lincs_drug_activity$standard_value < 100000 & lincs_drug_activity$pert_iname %in% lincs_drug_prediction_subset$pert_iname]

col_sorted = sort(cor(drug_dz_signature_rank, method="spearman")["log2FoldChange",-1])    
drug_dz_signature_rank = drug_dz_signature_rank[,c("log2FoldChange", names(col_sorted))]

drug_names = sapply(2:ncol(drug_dz_signature_rank), function(id){
  lincs_drug_prediction$pert_iname[paste("X",lincs_drug_prediction$id, sep="") == names(drug_dz_signature_rank)[id]][1]
})

#selected_drug_names

#drug_dz_signature_rank = drug_dz_signature_rank[, drug_names %in% c("log2FoldChange", selected_drug_names)]

pdf(paste( cancer, "/lincs_reverse_expression_selected.pdf", sep=""))
  colPal <- col
  par(mar=c(13, 6, 2, 0.5))
  axiscolor = sapply(c(cancer, as.character(drug_names)), function(name){
    if (name == cancer){
      "black"
    }else if (name %in% selected_drug_names){
      "red"
    }else{
      "green"
    }
  })
  image(t(drug_dz_signature_rank), col=colPal,   axes=F, srt=45)
  axis(1,  at= seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), labels= FALSE)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), c(-0.05),
       labels = c( cancer,as.character(drug_names)),col=axiscolor, srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=1.2)
dev.off()

