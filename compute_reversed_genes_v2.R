#find reversed genes by other anti-cancer drugs
#statistical test between active and inactive drugs
#leave one-drug-out cross valdiation
#revised: replace DER based on reviewers' comments. 
#use absolute position to replace DER
#include Achilles shRNA data to validate targets. This part was not added in the manuscript.

library(pheatmap)
library("gplots")
library(ConsensusClusterPlus)
library(metap)
library(RColorBrewer)
library(ggplot2)
require(Heatplus)

heatmap_col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)


cancer <- "LIHC"
landmark <- 1
fdr_cutoff <- 0.25
active_cutoff <- 10000
inactive_cutoff <- 10000
LIHC_validated_drugs <-  c("CGK-733", "strophanthidin", "FCCP",  "pyrvinium-pamoate")
LIHC_validated_drug_IC50 <-  c(3.18, 0.72, 1.78, 0.07)

#CMAP score output
output_path <- paste(cancer, "/lincs_score_", landmark, ".csv", sep="")
lincs_drug_prediction <- read.csv(output_path, stringsAsFactors = F)

if (cancer == "LIHC"){
  output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
  lincs_drug_prediction_all <- read.csv(output_path, stringsAsFactors = F)
  lincs_drug_prediction <- subset(lincs_drug_prediction_all, pert_iname %in% c(as.character(lincs_drug_prediction$pert_iname), as.character(LIHC_validated_drugs)))
  lincs_drug_prediction$RGES <- lincs_drug_prediction$cmap_score
}


load("raw/lincs/lincs_signatures_cmpd_landmark.RData")


########
load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #
res$GeneID <- sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol <- sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature <- subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf )
dz_sig <- subset(dz_signature, select=c("GeneID", "log2FoldChange"))

##########
#visualized gene reversed
#only pick the signatures from close to the median
cell_lines <- read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))
#choose profiles from the same lineage
lincs_drug_prediction_subset <- subset(lincs_drug_prediction, cell_id %in% cell_lines$LINCS, select=c("id", "RGES", "pert_iname", "pert_dose", "pert_time"))

###selecting median still sounds odd.
drug_cmap_score <- aggregate(RGES ~ pert_iname, lincs_drug_prediction_subset, median)
drug_instances_median  <- merge(lincs_drug_prediction_subset, drug_cmap_score, by = c("pert_iname"))
drug_instances_median$diff <- abs(drug_instances_median$RGES.x - drug_instances_median$RGES.y) #cmap_score.y is the median
drug_instances_min_diff <- aggregate(diff ~ pert_iname, drug_instances_median, min)
drug_instances_select <- merge(drug_instances_median, drug_instances_min_diff, by=c("pert_iname", "diff"))
drug_instances_select <- drug_instances_select[!duplicated(drug_instances_select$pert_iname), ]
sig_id_selects <- drug_instances_select$id

drug_dz_signature <- merge(dz_sig, data.frame(GeneID = rownames(lincs_signatures), lincs_signatures[, as.character(sig_id_selects)]),  by="GeneID", suffixes='')

#########################
###
#visualize the reversed gene expression
#reorder drug dz signatures
gene_ids <- drug_dz_signature$GeneID
drug_dz_signature_rank <- drug_dz_signature[,-1]
for (i in 1:ncol(drug_dz_signature_rank)){
  drug_dz_signature_rank[,i] <- rank(-1 * drug_dz_signature_rank[,i] ) #highly expressed genes ranked on the top
}
gene_ids_rank <- gene_ids[order(drug_dz_signature_rank[,1])]
drug_dz_signature_rank <- drug_dz_signature_rank[order(drug_dz_signature_rank[,1]),] #order by disease expression

lincs_drug_activity <- read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity <- unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",  	"cell_line")))
lincs_drug_activity <- aggregate(standard_value ~ pert_iname, lincs_drug_activity, median)
#active_cutoff <- inactive_cutoff <- median(lincs_drug_activity$standard_value)
active_drug_names <- lincs_drug_activity$pert_iname[lincs_drug_activity$standard_value <= active_cutoff & lincs_drug_activity$pert_iname %in% lincs_drug_prediction_subset$pert_iname]
inactive_drug_names <- lincs_drug_activity$pert_iname[lincs_drug_activity$standard_value > inactive_cutoff & lincs_drug_activity$pert_iname %in% lincs_drug_prediction_subset$pert_iname]

if (cancer == "LIHC"){active_drug_names <- c(active_drug_names, LIHC_validated_drugs)}

col_sorted <- sort(cor(drug_dz_signature_rank, method="spearman")["log2FoldChange",-1])    
drug_dz_signature_rank <- drug_dz_signature_rank[,c("log2FoldChange", names(col_sorted))]

drug_names <- sapply(2:ncol(drug_dz_signature_rank), function(id){
  lincs_drug_prediction$pert_iname[paste("X",lincs_drug_prediction$id, sep="") == names(drug_dz_signature_rank)[id]][1]
})



#####################
#find reversed genes using DER concept
###
#leave one drug out
#
###
#select reversed genes based on enrichment 
up_genes <- dz_sig$GeneID[dz_sig$GeneID %in% rownames(lincs_signatures) & dz_sig$log2FoldChange > 0]
down_genes <- dz_sig$GeneID[dz_sig$GeneID %in% rownames(lincs_signatures) & dz_sig$log2FoldChange < 0]

all_drugs <- unique(c(active_drug_names, inactive_drug_names))
down_p_value_matrix <- matrix(NA, nrow=length(all_drugs), ncol=length(up_genes))
up_p_value_matrix <- matrix(NA, nrow=length(all_drugs), ncol=length(down_genes))

drug_signature_rank <- lincs_signatures[, as.character(sig_id_selects)]
for (i in 1:ncol(drug_signature_rank)){
  drug_signature_rank[,i] <- rank(-1 * drug_signature_rank[,i] ) #highly expressed genes ranked on the top
}

for (i in 1:length(all_drugs)){
  #print(i)
  drug_left <- all_drugs[i]
  
  effective_drug_signature_ids <- lincs_drug_prediction[lincs_drug_prediction$pert_iname %in% active_drug_names & !lincs_drug_prediction$pert_iname %in% drug_left, "id"]
  ineffective_drug_signature_ids <- lincs_drug_prediction[lincs_drug_prediction$pert_iname %in% inactive_drug_names & !lincs_drug_prediction$pert_iname %in% drug_left, "id"]
  
  up_ps <- NULL
  for (gene in up_genes){
    effective_ranks <- as.numeric(drug_signature_rank[gene, colnames(drug_signature_rank) %in% effective_drug_signature_ids])
    ineffective_ranks <- as.numeric(drug_signature_rank[gene, colnames(drug_signature_rank) %in% ineffective_drug_signature_ids])
    up_ps <- c(up_ps, wilcox.test(effective_ranks, ineffective_ranks, alternative = "greater")$p.value)
  }
  p_adj <- p.adjust(up_ps, method="fdr")
  down_p_value_matrix[i,] <- p_adj
  
  down_ps <- NULL
  for (gene in down_genes){
    effective_ranks <- as.numeric(drug_signature_rank[gene, colnames(drug_signature_rank) %in% effective_drug_signature_ids])
    ineffective_ranks <- as.numeric(drug_signature_rank[gene, colnames(drug_signature_rank) %in% ineffective_drug_signature_ids])
    down_ps <- c(down_ps, wilcox.test(effective_ranks, ineffective_ranks, alternative = "less")$p.value)
  }
  p_adj <- p.adjust(down_ps, method="fdr")
  up_p_value_matrix[i,] <- p_adj
  
}

suppress_gene_p <- data.frame(GeneID = up_genes, sig_down_count = apply(down_p_value_matrix, 2, function(values){sum(values<fdr_cutoff)}),
                    sig_down_meta = apply(down_p_value_matrix, 2, function(values){sumlog(values)$p})
                    )
suppress_gene_p_annot <- merge(suppress_gene_p, res, by.x="GeneID", by.y="GeneID")
suppress_gene_p_annot <- suppress_gene_p_annot[order(suppress_gene_p_annot$sig_down_meta), ]
suppress_gene_p_annot$sig_down_meta_padj <- p.adjust(suppress_gene_p_annot$sig_down_meta, method = "fdr")

induce_gene_p <- data.frame(GeneID = down_genes, sig_up_count = apply(up_p_value_matrix, 2, function(values){sum(values<fdr_cutoff)}),
                             sig_up_meta = apply(up_p_value_matrix, 2, function(values){sumlog(values)$p})
)
induce_gene_p_annot <- merge(induce_gene_p, res, by.x="GeneID", by.y="GeneID")
induce_gene_p_annot <- induce_gene_p_annot[order(induce_gene_p_annot$sig_up_meta), ]
induce_gene_p_annot$sig_up_meta_padj <- p.adjust(induce_gene_p_annot$sig_up_meta, method = "fdr")

write.table(induce_gene_p_annot, paste(cancer,"/target_reversed_score_induced.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)
write.table(suppress_gene_p_annot, paste(cancer,"/target_reversed_score_suppressed.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)


#visualize only reversed geneså
###genes reversed
drug_signature_rank_ratio <- drug_signature_rank[c(as.character(induce_gene_p_annot$GeneID[induce_gene_p_annot$sig_up_count == length(all_drugs)]),
                                                  as.character(suppress_gene_p_annot$GeneID[suppress_gene_p_annot$sig_down_count == length(all_drugs)])), ] 
drug_signature_rank_ratio <- drug_signature_rank_ratio / 978

gene_ids <- rownames(drug_signature_rank_ratio)
gene_ids_annot <- merge(gene_ids, res, by.x=1, by.y="GeneID", sort=F)

drug_names <- sapply(1:ncol(drug_signature_rank_ratio), function(id){
  lincs_drug_prediction$pert_iname[paste(lincs_drug_prediction$id, sep="") == as.character(colnames(drug_signature_rank_ratio)[id])][1]
})

rownames(drug_signature_rank_ratio) <- gene_ids_annot$symbol
colnames(drug_signature_rank_ratio) <- drug_names

annotation <- rbind(data.frame(drug = active_drug_names, type = "effective"), data.frame(drug = inactive_drug_names, type = "ineffective"))
rownames(annotation) <- annotation[,1]
annotation$type <- factor(annotation$type)
annotation <- subset(annotation, select="type")
names(annotation) <- "Activity"
Var1 <- c("green", "red")
names(Var1) <- c("effective", "ineffective")
ann_colors <- list(Activity = c(effective = "green", ineffective="blue"))

my.cols <- heatmap_col # redgreen(100) # brewer.pal(9, "Blues")
pheatmap((drug_signature_rank_ratio), annotation = annotation,  annotation_legend = F, col = my.cols, 
         cellheight = 14, cellwidth = 1, show_rownames = T, legend=F,  show_colnames=F,
         clustering_distance_cols = "euclidean", annotation_colors = ann_colors,
         filename = paste( "fig/targets_", cancer, "_reversed_genes_v2.pdf", sep="")) #

#order by activity
rownames(lincs_drug_activity) <- lincs_drug_activity$pert_iname
if (cancer == "LIHC"){
  LIHC_validated <- data.frame(pert_iname = LIHC_validated_drugs, standard_value = LIHC_validated_drug_IC50 * 1000)
  rownames(LIHC_validated) <- LIHC_validated$pert_iname
  lincs_drug_activity <- rbind(lincs_drug_activity, LIHC_validated)
}
drug_activity <- lincs_drug_activity[colnames(drug_signature_rank_ratio), ]
drug_signature_rank_ratio <- drug_signature_rank_ratio[, drug_activity$pert_iname[order(drug_activity$standard_value)]]
pheatmap((drug_signature_rank_ratio), annotation = annotation,  annotation_legend = F, col = my.cols, 
         cellheight = 14, cellwidth = 1, show_rownames = T, legend=F,  show_colnames=F, cluster_cols=F, treeheight_row = 1, 
         clustering_distance_cols = "euclidean", annotation_colors = ann_colors, 
         filename = paste( "fig/targets_", cancer, "_reversed_genes_order_by_activity_v2.pdf", sep="")) #

pheatmap((drug_signature_rank_ratio), annotation = annotation,  annotation_legend = F, col = my.cols, 
         cellheight = 14, cellwidth = 1, show_rownames = T, legend=T,  show_colnames=F, cluster_cols=F, treeheight_row = 1, 
         clustering_distance_cols = "euclidean", annotation_colors = ann_colors, 
         filename = paste( "fig/targets_", cancer, "_reversed_genes_order_by_activity_with_legend_v2.pdf", sep="")) #


col_annot <- data.frame(IC50 = log(lincs_drug_activity[colnames(drug_signature_rank_ratio), "standard_value"], 10))
rownames(col_annot) <- colnames(drug_signature_rank_ratio)

'reg1 = annHeatmap2(drug_signature_rank_ratio,                   
                   ann=list(Col=list(data=col_annot[colnames(drug_signature_rank_ratio),])), 
                   scale="none",
                   cluster=list(Col=list(cuth=0.3)),
                   labels=list(Col = list(nrow=1,cex=0.01)),
                   legend=1)

plot(reg1)'

out <- pheatmap((drug_signature_rank_ratio), annotation = annotation,  annotation_legend = T, col = my.cols, 
                cellheight = 12, cellwidth = 8, show_rownames = T, legend=T,  
                clustering_distance_cols = "euclidean", annotation_colors = ann_colors)

pheatmap_mat <- drug_signature_rank_ratio[c(out$tree_row[["order"]]),out$tree_col[["order"]]]

#visualize significant genes by fold change

pdf(paste("fig/fc_", cancer, ".pdf", sep=""))
res_reversed <- res[res$symbol %in% rownames(pheatmap_mat), ]
rownames(res_reversed) <- res_reversed$symbol
res_reversed <- res_reversed[rownames(pheatmap_mat), ]
res_reversed$symbol <- factor(res_reversed$symbol, levels = rev(res_reversed$symbol))
res_reversed$up_down <- ifelse(res_reversed$log2FoldChange>0, "up", "down")

ggplot(res_reversed, aes(x=(symbol), y = log2FoldChange, fill = up_down)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), axis.text = element_text(colour = "black", size=10)) +
  xlab("") +
  ylab("") +  #log 2 fold change (tumor vs. non tumor)
  scale_fill_manual(values=c(  "green","red")) +
  theme(axis.title.x = element_text(size=17),axis.text.x = element_text(size=17), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y=element_blank(),axis.line.y = element_blank(),
        plot.margin = unit(c(8, 0.5, 0.5, 13), "cm"),
        legend.position="none")  ##top, right, bottom, left
dev.off()



####
#valdiate shRNA data
###
shrna <- read.table("~/Documents/stanford/tumor_cell_line/data/Achilles_QC_v2.4.3.rnai.Gs.gct", sep="\t", skip=2, header=T)
#shrna <- read.table("~/Downloads/Achilles_v3.3.8.Gs.gct", sep="\t", skip=2, header=T)

#shrna <- shrna[!duplicated(shrna$Description),]
#rownames(shrna) <- shrna$Description

if (cancer == "LIHC"){
  tissue <- "LIVER"
}else if (cancer == "COAD"){
  tissue <- "LARGE_INTESTINE"
}else if (cancer == "BRCA" | cancer == "ER"){
  tissue <- "BREAST"
}
cell_lines <- colnames(shrna)[grep(tissue, colnames(shrna))] #HT29; CAL120_BREAST T47D_BREAST LARGE_INTESTINE

all_viability <- shrna[!shrna$Description %in% suppress_gene_p_annot$symbol, cell_lines]
suppress_shrna_viability <- shrna[shrna$Description %in% suppress_gene_p_annot$symbol[suppress_gene_p_annot$sig_down_count == length(all_drugs)], cell_lines]
shrna_viability <- shrna[shrna$Description %in% suppress_gene_p_annot$symbol[suppress_gene_p_annot$sig_down_count < length(all_drugs)], cell_lines]
induce_shrna_viability <- shrna[shrna$Description %in% induce_gene_p_annot$symbol[induce_gene_p_annot$sig_up_count == length(all_drugs)], cell_lines]
mean(as.matrix(induce_shrna_viability))
mean(as.matrix(suppress_shrna_viability))
mean(as.matrix(shrna_viability))

t.test(as.numeric(as.matrix(suppress_shrna_viability)), as.numeric(as.matrix(shrna_viability)))
t.test(as.numeric(as.matrix(suppress_shrna_viability)), as.numeric(as.matrix(all_viability)))

boxplot(as.numeric(as.matrix(suppress_shrna_viability)), 
        as.numeric(as.matrix(shrna_viability)),
        as.numeric(as.matrix(all_viability)))

#####################
#enrichment of a single genes
gene <- "4172"
colPal <- redgreen(100)

drug_signature_enriched <- matrix(0, nrow = ncol(drug_signature_rank), ncol = 978 )
drug_names <- sapply(1:ncol(drug_signature_rank), function(id){
  lincs_drug_prediction$pert_iname[paste(lincs_drug_prediction$id, sep="") == colnames(drug_signature_rank)[id]][1]
})
colnames(drug_signature_rank) <- drug_names
row.names(drug_signature_enriched) <- drug_names
gene_rank <- drug_signature_rank[as.character(gene), ]
drug_signature_enriched <- drug_signature_enriched[names(sort(gene_rank)), ]
drug_signature_enriched <- rbind(drug_signature_enriched[rownames(drug_signature_enriched) %in% inactive_drug_names, ], 
                                drug_signature_enriched[rownames(drug_signature_enriched) %in% active_drug_names,])

for (i in 1:nrow(drug_signature_enriched)){
  drug_signature_enriched[i, round(gene_rank[rownames(drug_signature_enriched)[i]])] <- 1
}

pdf(paste(cancer, "/reverse_", gene, ".pdf", sep=""))
par(mar=c(0.5, 10, 0.5, 0.5))
image(t(drug_signature_enriched),  axes=F, srt=45, col=heat.colors(12))
axis(2,  at=seq(0,1,length.out= nrow( drug_signature_enriched ) ), labels= F)
text(y = seq(0,1,length.out = nrow( drug_signature_enriched )  ), c(-0.05),
     col = c(rep("blue", length(inactive_drug_names)), rep("green", length(active_drug_names))), 
     labels = rownames( drug_signature_enriched ) , srt = 0, pos=2, offset=0.05, xpd = TRUE, cex=1.2)
dev.off()

drug_signature_enriched_merged <- rbind(apply(drug_signature_enriched[inactive_drug_names, ], 2, sum),
                                       apply(drug_signature_enriched[active_drug_names, ], 2, sum))
rownames(drug_signature_enriched_merged) <- c("ineffective", "effective")

drug_signature_enriched_merged[drug_signature_enriched_merged>0] <- 1
wilcox.test(which(drug_signature_enriched_merged["effective",] == 1),
            which(drug_signature_enriched_merged["ineffective",] == 1))

#colar validated drugs
LIHC_drug_pos <- c(which(drug_signature_enriched["CGK-733",] == 1), which(drug_signature_enriched["strophanthidin",] == 1), 
                  which(drug_signature_enriched["FCCP",] == 1), which(drug_signature_enriched["pyrvinium-pamoate",] == 1))

drug_signature_enriched_merged["effective", LIHC_drug_pos] <- 2

pdf(paste(cancer, "/reverse_", gene, "_merged.pdf", sep=""))
par(mar=c(30, 10, 0.5, 0.5))
image(t(drug_signature_enriched_merged),  col=  colorpanel(3, low="gray", mid="black", high = "red"), axes=F, srt=45)
#text(y = seq(0,1,length.out = nrow( drug_signature_enriched_merged )  ), c(-0.05),
#     labels= NULL, srt = 0, pos=2, offset=0.05, xpd = TRUE, cex=1.5) #labels = rownames( drug_signature_enriched_merged )
dev.off()


##########
#unsupervise clustering to find common clusters
#drugs with distinct mechanism may be located in different clusters
drug_names <- sapply(1:ncol(drug_signature_rank), function(id){
  lincs_drug_prediction$pert_iname[paste(lincs_drug_prediction$id, sep="") == as.character(colnames(drug_signature_rank)[id])][1]
})
colnames(drug_signature_rank) <- drug_names
drug_signature_rank_effective <- drug_signature_rank[, colnames(drug_signature_rank) %in% c(active_drug_names, inactive_drug_names) ]

results <- ConsensusClusterPlus(as.matrix(drug_signature_rank_effective),maxK=6,reps=50,pItem=0.8,pFeature=1,
                               title=paste(cancer, "_cluster",sep=""),clusterAlg="hc",distance="euclidean",seed=1262118388.71279,plot="png")
sort(results[[2]][["consensusClass"]])

reversed_genes <- res$symbol[res$GeneID %in% c(up_genes, down_genes)]

clust_member <- results[[2]][["consensusClass"]]
sum(names(clust_member[clust_member == 1]) %in% inactive_drug_names)
sum(names(clust_member[clust_member == 2]) %in% active_drug_names)

#####################
#retrieve mutation profiles
#deprecated; using data from cbioportal and visualize using oncoprint.
####################
reversed_genes <- rownames(res_reversed)
library(cgdsr)
mycgds <- CGDS("http://www.cbioportal.org/public-portal/")
# Get list of cancer studies at server
a <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudies <- getCancerStudies(mycgds)[,1]

if (cancer == "COAD"){
  mycancerstudy <- "coadread_tcga_pub"
}else if (cancer == "LIHC"){
  mycancerstudy <- "lihc_tcga"
}else if (cancer == "BRCA"){
  mycancerstudy <- "brca_tcga"
}
# Get available genetic profiles
mycaselist <- getCaseLists(mycgds,mycancerstudy)[1,1]
mygeneticprofiles <- getGeneticProfiles(mycgds,mycancerstudy)[,1]
######################

####
#mutation profiles
mygeneticprofile <- mygeneticprofiles[grep("mutations", mygeneticprofiles)]
mutations <- getProfileData(mycgds,reversed_genes ,mygeneticprofile,mycaselist)
mutations_ratio <- apply(mutations, 2, function(x){sum(!is.na(x) & x != "NaN")/length(x)})

mutations_ratio <- data.frame(mutations_ratio)
mutations_ratio$symbol <- rownames(mutations_ratio)

ggplot(mutations_ratio, aes(x=(symbol), y = mutations_ratio)) +
  geom_bar(stat = "identity", fill="green") + 
  coord_flip() +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), axis.text = element_text(colour = "black", size=10)) +
  xlab("") + ylim(0, 0.4) + 
  ylab("mutation") + 
  scale_fill_manual(values=c( "red", "green")) +
  theme(axis.text.x = element_text(size=17), axis.text.y = element_text(size=17), axis.title.y = element_text(size=18),  axis.title.x = element_text(size=18)) 

#copy number variation
mygeneticprofile <- mygeneticprofiles[grep("gistic", mygeneticprofiles)]

#correct symbols
reversed_genes[reversed_genes == "GPER"] <- "GPER1"

cnv <- getProfileData(mycgds,reversed_genes ,mygeneticprofile,mycaselist)
cnv <- cnv[, reversed_genes]
cnv_deletion <- apply(cnv, 2, function(x){sum(x < -1, na.rm = T)/length(x)})
cnv_amplication <- apply(cnv, 2, function(x){sum(x > 1, na.rm = T)/length(x)})

cnv_ratio <- cbind(data.frame(cnv_deletion), data.frame(cnv_amplication))
cnv_ratio$symbol <- rownames(cnv_ratio)
cnv_ratio$del_amp <- "Deletion"
cnv_ratio$del_amp[cnv_ratio$cnv_deletion < cnv_ratio$cnv_amplication] <- "Amplification"
cnv_ratio$ratio <- sapply(1:nrow(cnv_ratio), function(x){max(cnv_ratio$cnv_deletion[x], cnv_ratio$cnv_amplication[x])})

cnv_ratio$symbol <- factor(cnv_ratio$symbol, levels = rev(cnv_ratio$symbol))
pdf(paste("fig/cn_", cancer, ".pdf", sep=""))
ggplot(cnv_ratio, aes(x=(symbol), y = ratio, fill = del_amp)) +
  geom_bar(stat = "identity") + 
  coord_flip() +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), axis.text = element_text(colour = "black", size=10)) +
  xlab("") + ylim(0, 0.3) + 
  ylab("Amplification/Deletion") + 
  scale_fill_manual(values=c( "red", "green")) +
  theme(axis.title.x = element_text(size=17),axis.text.x = element_text(size=17), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y=element_blank(),axis.line.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position="none")  ##top, right, bottom, left
dev.off()


###################
###using SAM to compute differential expression between effective group and ineffective group.
##################
drug_signature <- lincs_signatures[, as.character(sig_id_selects)]

drug_names_no_ranked <- sapply(1:ncol(drug_signature), function(id){
  lincs_drug_prediction$pert_iname[paste(lincs_drug_prediction$id, sep="") == colnames(drug_signature)[id]][1]
})

comparison_frame <- drug_signature
gene_ids <- rownames(drug_signature)
sample_class <- rep(0, length(drug_names_no_ranked))
sample_class[drug_names_no_ranked %in% active_drug_names] <- 1
method <- "SAM"
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


#visualize
comparison_frame_subset <- comparison_frame
rownames(comparison_frame_subset) <- gene_ids
comparison_frame_subset <- comparison_frame_subset[rownames(comparison_frame_subset) %in% c(genes.up$probe[genes.up$q.value < 0.25], genes.down$probe[genes.down$q.value < 0.25]), ]

drug_names <- sapply(1:ncol(comparison_frame_subset), function(id){
  lincs_drug_prediction$pert_iname[paste(lincs_drug_prediction$id, sep="") == as.character(colnames(comparison_frame_subset)[id])][1]
})

rownames(comparison_frame_subset) <- merge(data.frame(GeneID=rownames(comparison_frame_subset)), res, by = "GeneID", sort=F)$symbol
colnames(comparison_frame_subset) <- drug_names

annotation <- rbind(data.frame(drug = active_drug_names, type = "effective"), data.frame(drug = inactive_drug_names, type = "ineffective"))
rownames(annotation) <- annotation[,1]
annotation$type <- factor(annotation$type)
annotation <- subset(annotation, select="type")
names(annotation) <- "Activity"
Var1 <- c("green", "red")
names(Var1) <- c("effective", "ineffective")
ann_colors <- list(Activity = c(effective = "green", ineffective="blue"))

my.cols <- greenred(100) # brewer.pal(9, "Blues")
pheatmap((comparison_frame_subset), annotation = annotation,  annotation_legend = T, col = my.cols, 
         cellheight = 12, cellwidth = 8, show_rownames = T, legend=T,  
         clustering_distance_cols = "euclidean", annotation_colors = ann_colors, 
         filename = paste( "fig/targets_", cancer, "_reversed_genes_SAM.pdf", sep="")) #

######
#correlate rank with efficacy
######
drug_signature <- lincs_signatures[, as.character(sig_id_selects)]
for (i in 1:ncol(drug_signature_rank)){
  drug_signature_rank[,i] <- rank(-1 * drug_signature[,i] ) #highly expressed genes ranked on the top
}
drug_names <- sapply(1:ncol(drug_signature_rank), function(id){
  lincs_drug_prediction$pert_iname[paste(lincs_drug_prediction$id, sep="") == colnames(drug_signature_rank)[id]][1]
})
rownames(lincs_drug_activity) <- lincs_drug_activity$pert_iname

cor_rank_activity <- cor(lincs_drug_activity[drug_names, "standard_value"], t(drug_signature), method = "spearman")

cor_rank_activity_res <- merge(data.frame(GeneID = colnames(cor_rank_activity), cor = as.numeric(cor_rank_activity)), res, by = "GeneID")
cor_rank_activity_res <- cor_rank_activity_res[order(cor_rank_activity_res$cor), ]

plot(drug_signature["4582",], log(lincs_drug_activity[drug_names, "standard_value"], 10))



wilcox.test(as.numeric(as.matrix(suppress_shrna_viability)), as.numeric(as.matrix(all_viability)))

shrna[shrna$Description == "EGR1", cell_lines]


