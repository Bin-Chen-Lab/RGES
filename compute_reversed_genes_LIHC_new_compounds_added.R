#find reversed genes by other anti-cancer drugs
#t test between active and inactive drugs
#leave one-drug-out cross valdiation
library(pheatmap)
library("gplots")
library(ConsensusClusterPlus)
library(metap)
library(RColorBrewer)


cancer = "ER"
landmark = 1
fdr_cutoff = 0.25

validated_drugs =  c("CGK-733", "strophanthidin", "FCCP",  "pyrvinium-pamoate")
#CMAP score output
output_path <- paste(cancer, "/lincs_score_", landmark, ".csv", sep="")
lincs_drug_prediction_chembl = read.csv(output_path, stringsAsFactors = F)

output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
lincs_drug_prediction = read.csv(output_path, stringsAsFactors = F)
lincs_drug_prediction = subset(lincs_drug_prediction, pert_iname %in% c(lincs_drug_prediction_chembl$pert_iname, validated_drugs))
lincs_drug_prediction$RGES = lincs_drug_prediction$cmap_score

load("raw/lincs/lincs_signatures_cmpd_landmark.RData")

########
load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #
res$GeneID = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf )
dz_sig = subset(dz_signature, select=c("GeneID", "log2FoldChange"))

##########
#visualized gene reversed
#only pick the signatures from close to the median
cell_lines = read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))
#choose profiles from the same lineage
lincs_drug_prediction_subset = subset(lincs_drug_prediction, cell_id %in% cell_lines$LINCS, select=c("id", "RGES", "pert_iname", "pert_dose", "pert_time"))

###selecting median still sounds odd.
drug_cmap_score = aggregate(RGES ~ pert_iname, lincs_drug_prediction_subset, median)
drug_instances_median  = merge(lincs_drug_prediction_subset, drug_cmap_score, by = c("pert_iname"))
drug_instances_median$diff = abs(drug_instances_median$RGES.x - drug_instances_median$RGES.y) #cmap_score.y is the median
drug_instances_min_diff = aggregate(diff ~ pert_iname, drug_instances_median, min)
drug_instances_select = merge(drug_instances_median, drug_instances_min_diff, by=c("pert_iname", "diff"))
drug_instances_select = drug_instances_select[!duplicated(drug_instances_select$pert_iname), ]
sig_id_selects = drug_instances_select$id

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

active_drug_names = c(active_drug_names, validated_drugs)

col_sorted = sort(cor(drug_dz_signature_rank, method="spearman")["log2FoldChange",-1])    
drug_dz_signature_rank = drug_dz_signature_rank[,c("log2FoldChange", names(col_sorted))]

drug_names = sapply(2:ncol(drug_dz_signature_rank), function(id){
  lincs_drug_prediction$pert_iname[paste("X",lincs_drug_prediction$id, sep="") == names(drug_dz_signature_rank)[id]][1]
})

pdf(paste( "fig/", cancer, "_lincs_reverse_expression.pdf", sep=""))
  colPal <- greenred(100)
  par(mar=c(13, 6, 2, 0.5))
  axiscolor = sapply(c(cancer, as.character(drug_names)), function(name){
    if (name == cancer){
      "black"
    }else if (name %in% active_drug_names){
      "red"
    }else{
      "green"
    }
  })
  image(t(drug_dz_signature_rank), col=colPal,   axes=F, srt=45)
  axis(1,  at= seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), labels= FALSE)
  text(x = seq(0,1,length.out=ncol( drug_dz_signature_rank ) ), c(-0.05),
       labels = c( cancer,as.character(drug_names)),col=axiscolor, srt = 45, pos=2,offset=0.05, xpd = TRUE, cex=0.4)
dev.off()

##########
#unsupervise clustering to find common clusters
# results = ConsensusClusterPlus(as.matrix(drug_dz_signature_reversed),maxK=6,reps=50,pItem=0.8,pFeature=1,
#                                title=paste(cancer, "_cluster",sep=""),clusterAlg="hc",distance="euclidean",seed=1262118388.71279,plot="png")
# sort(results[[2]][["consensusClass"]])

#####################
#find reversed genes
###
#leave one drug out
all_drugs = unique(active_drug_names, inactive_drug_names)
down_p_value_matrix = matrix(NA, nrow=length(all_drugs), ncol=nrow(drug_dz_signature_rank))
up_p_value_matrix = matrix(NA, nrow=length(all_drugs), ncol=nrow(drug_dz_signature_rank))

for (i in 1:length(all_drugs)){
  print(i)
  drug_left = all_drugs[i]
  active_drug_names_new = active_drug_names[!active_drug_names %in% drug_left]
  inactive_drug_names_new = inactive_drug_names[!inactive_drug_names %in% drug_left]
  
  lincs_sig_info_id = paste("X", lincs_drug_prediction$id, sep="")
  active_drugs = lincs_sig_info_id[lincs_drug_prediction$pert_iname %in% active_drug_names_new & lincs_sig_info_id %in% colnames(drug_dz_signature_rank)]
  inactive_drugs = lincs_sig_info_id[lincs_drug_prediction$pert_iname %in% inactive_drug_names_new & lincs_sig_info_id %in% colnames(drug_dz_signature_rank)]

  down_p_value = sapply(1:nrow(drug_dz_signature_rank), function(id){
    t.test(as.numeric(drug_dz_signature_rank[id, active_drugs] - drug_dz_signature_rank[id,1]) , as.numeric(drug_dz_signature_rank[id, inactive_drugs] - drug_dz_signature_rank[id,1]), alternative=c("greater"))$p.value
  })
  up_p_value = sapply(1:nrow(drug_dz_signature_rank), function(id){
    t.test(as.numeric(drug_dz_signature_rank[id, active_drugs] - drug_dz_signature_rank[id,1]) , as.numeric(drug_dz_signature_rank[id, inactive_drugs] - drug_dz_signature_rank[id,1]), alternative=c("less"))$p.value
  })
  
  down_p_value_matrix[i,] =  p.adjust(down_p_value, method="fdr") #down_p_value #
  up_p_value_matrix[i,] =  p.adjust(up_p_value, method="fdr") #up_p_value #
}

gene_p = data.frame(gene = gene_ids_rank, sig_down_count = apply(down_p_value_matrix, 2, function(values){sum(values<0.25)}),
                    sig_up_count = apply(up_p_value_matrix, 2, function(values){sum(values<0.25)}),
                    sig_down_meta = apply(down_p_value_matrix, 2, function(values){sumlog(values)$p}),
                    sig_up_meta = apply(up_p_value_matrix, 2, function(values){sumlog(values)$p})
                    )
gene_p_annot = merge(gene_p, res, by.x="gene", by.y="GeneID")
gene_p_annot = gene_p_annot[order(gene_p_annot$sig_down_meta), ]
gene_p_annot$sig_up_meta_padj = p.adjust(gene_p_annot$sig_up_meta, method="fdr")
gene_p_annot$sig_down_meta_padj = p.adjust(gene_p_annot$sig_down_meta, method = "fdr")

#valid_genes = data.frame(gene=gene_ids_rank, up_p=up_p_value, up_q= up_fdr, down_p = down_p_value, down_q = down_fdr)
#valid_genes = valid_genes[order(valid_genes$down_q),]
#valid_genes = subset(valid_genes, down_q < fdr_cutoff | up_q < fdr_cutoff)
#valid_genes = merge(valid_genes, res, by.x="gene", by.y="GeneID")

#write.table(valid_genes[valid_genes$up_q < fdr_cutoff,], paste(cancer,"/target_up_genes_classify_good_bad_drug.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
#write.table(valid_genes[valid_genes$down_q < fdr_cutoff,], paste(cancer,"/target_down_genes_classify_good_bad_drug.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
write.table(gene_p_annot, paste(cancer,"/target_reversed_score.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)


#visualize only reversed geneså
###genes reversed
gene_ids = drug_dz_signature$GeneID
gene_ids_annot <- merge(gene_ids, res, by.x=1, by.y="GeneID", sort=F)

drug_dz_signature_rank = drug_dz_signature[,-1]
for (i in 1:ncol(drug_dz_signature_rank)){
  drug_dz_signature_rank[,i] = rank(1 * drug_dz_signature_rank[,i] ) #highly expressed genes ranked on the top
}

drug_names = sapply(2:ncol(drug_dz_signature_rank), function(id){
  lincs_drug_prediction$pert_iname[paste("X",lincs_drug_prediction$id, sep="") == as.character(names(drug_dz_signature_rank)[id])][1]
})

drug_dz_signature_reversed = (drug_dz_signature_rank - drug_dz_signature_rank[, 1])
drug_dz_signature_reversed = drug_dz_signature_reversed[, -c(1)]
rownames(drug_dz_signature_reversed) = gene_ids_annot$symbol
colnames(drug_dz_signature_reversed) = drug_names
drug_dz_signature_reversed = drug_dz_signature_reversed[rownames(drug_dz_signature_reversed) %in% as.character(gene_p_annot$symbol[gene_p_annot$sig_down_count >= length(all_drugs) | gene_p_annot$sig_up_count >= length(all_drugs)]),]

annotation = rbind(data.frame(drug = active_drug_names, type = "effective"), data.frame(drug = inactive_drug_names, type = "ineffective"))
rownames(annotation) = annotation[,1]
annotation$type = factor(annotation$type)
annotation = subset(annotation, select="type")
names(annotation) = "Activity"
Var1 = c("green", "red")
names(Var1) = c("effective", "ineffective")
ann_colors = list(Activity = c(effective = "green", ineffective="blue"))

my.cols <- greenred(100) # brewer.pal(9, "Blues")
pheatmap((drug_dz_signature_reversed), annotation = annotation,  annotation_legend = T, col = my.cols, 
         cellheight = 12, cellwidth = 8, show_rownames = T, legend=T,  
         clustering_distance_cols = "euclidean", annotation_colors = ann_colors,
         filename = paste( "fig/targets_", cancer, "_reversed_genes.pdf", sep="")) #


##############
#compute using rank
drug_names_no_ranked = sapply(2:ncol(drug_dz_signature), function(id){
  lincs_drug_prediction$pert_iname[paste("X",lincs_drug_prediction$id, sep="") == names(drug_dz_signature_rank)[id]][1]
})

comparison_frame = drug_dz_signature[,-c(1)]
gene_ids = drug_dz_signature[, 1]
sample_class = rep(0, length(drug_names_no_ranked))
sample_class[drug_names_no_ranked %in% active_drug_names] = 1
method = "SAM"
if (method == 'RANKPROD_SINGLE') {
  library(RankProd)
  RP_result <- RP(comparison_frame, sample_class,gene.names=gene_ids, num.perm = 100, logged = F, na.rm = TRUE, plot = FALSE, rand = 123)
  # Leave logged=FALSE because topGene() converts the fold-change incorrectly!
  siggenes <- topGene(RP_result,cutoff=0.1,method="pfp",logged=FALSE,gene.names=gene_ids)
  # Normalize the results across methods
  siggenes.result = list()
  siggenes.result$UP <- (siggenes$Table1[,3:5])
  colnames(siggenes.result$UP) <- c("fold.change","q.value","p.value")
  siggenes.result$DOWN <- siggenes$Table2[,3:5]
  colnames(siggenes.result$DOWN) <- c("fold.change","q.value","p.value")
  # Since RANKPROD does the goofy condition 1 / condition 2, inverse the fold-change and convert to log if not logged. 
  siggenes.result$UP[,"fold.change"] <- log2(1/siggenes.result$UP[,"fold.change"])
  siggenes.result$DOWN[,"fold.change"] <- log2(1/siggenes.result$DOWN[,"fold.change"])
  
  siggenes.result$UP=data.frame(siggenes.result$UP,probe=rownames(siggenes.result$UP))
  siggenes.result$DOWN=data.frame(siggenes.result$DOWN,probe=rownames(siggenes.result$DOWN))
  
  
} else {
  library(siggenes)
  SAM_result <- sam(comparison_frame,sample_class,rand=123,gene.names=gene_ids) #q.version=1
  delta.table <- rbind(c(0,0,0),findDelta(SAM_result,fdr=0.9)) #fdr=0.9, too loose
  siggenes.table <- summary(SAM_result, delta.table[  dim(delta.table)[1],  1] );
  siggenes <- siggenes.table@mat.sig
  
  if( nrow(siggenes) > 0 ) {
    siggenes.result = list()
    siggenes.result$UP <- siggenes[siggenes$d.value > 0 & siggenes$rawp < 0.05,c("R.fold","q.value","rawp")]
    colnames(siggenes.result$UP) <- c("fold.change","q.value","p.value")
    siggenes.result$DOWN <- siggenes[siggenes$d.value < 0 & siggenes$rawp < 0.05,c("R.fold","q.value","rawp")]
    colnames(siggenes.result$DOWN) <- c("fold.change","q.value","p.value")
    
    
    siggenes.result$UP=data.frame(siggenes.result$UP,probe=rownames(siggenes.result$UP))
    siggenes.result$DOWN=data.frame(siggenes.result$DOWN,probe=rownames(siggenes.result$DOWN))
  }
  genes.up <- merge(siggenes.result$UP,res, by.x="probe",by.y="GeneID",sort=F,all.x=T)
  genes.down <- merge(siggenes.result$DOWN,res, by.x="probe",by.y="GeneID",sort=F,all.x=T)
}

write.table(rbind(genes.up,genes.down), paste(cancer,"/target_SAM_good_bad_drugs.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)


