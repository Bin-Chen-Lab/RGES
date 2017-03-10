setwd("~/Documents/stanford/tumor_cell_line/RGES_manuscript/data")

library(pheatmap)
library(gplots)
library("DESeq")
library("RColorBrewer")
my.cols <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256)) #RdYlBu


load("../../data/GSE62944/TumorCounts.RData")
load("../../data/GSE62944/NormalCounts.RData")
load("../../data/GSE62944/clinvar.RData")
load("../../data/GSE62944/NormalCancerType.RData")

cancer = "COAD"

load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC

#genes shared by landmark genes and dz signatures
landmark = read.csv("~/Desktop/lincs_landmark.csv") #get.landmark.info(con)

dz_padj_cutoff = 1E-3
dz_fc_cutoff = 1.5

dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < dz_padj_cutoff & abs(log2FoldChange) > dz_fc_cutoff & abs(log2FoldChange) != Inf )
dz_signature$symbol= sapply(dz_signature$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[1])
})
dz_signature$GeneID = sapply(dz_signature$id, function(id){
  as.character(unlist(strsplit(id, "\\|"))[2])
})

dz_signature = subset(dz_signature, GeneID %in% landmark$gene_id)


tumor_samples = table(clinvar$CancerType)
normal_samples = table(NormalCancerType$type)

TumorCounts_subset=TumorCounts[,rownames(clinvar[clinvar$CancerType == cancer,])]
NormalCounts_subset = NormalCounts[, NormalCancerType$sample[NormalCancerType$type == cancer]]
cds = newCountDataSet(cbind(TumorCounts_subset, NormalCounts_subset),  c(rep("tumor", ncol(TumorCounts_subset)), rep("non-tumor", ncol(NormalCounts_subset)))) #only integer is accepted
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds) #, method ="pooled", sharingMode="gene-est-only")
  
normalized_cds_all = counts( cds, normalized=TRUE )
normalized_cds_all = log(normalized_cds_all + 1)

#remove outlier samples
#normalized_cds_all = normalized_cds_all[, colnames(normalized_cds_all) %in% rownames(tumor_info[tumor_info$type == "non-tumor" | (tumor_info$type == "tumor" & tumor_info$quality != "outlier"),])]

normalized_cds_sig = normalized_cds_all[rownames(normalized_cds_all) %in% dz_signature$symbol,]

annotation = rbind(data.frame(tcga_barcode = colnames(TumorCounts_subset), type = "tumor"), data.frame(tcga_barcode = colnames(NormalCounts_subset), type = "normal") )
rownames(annotation) = annotation$tcga_barcode
annotation$type = as.factor(annotation$type)
annotation = subset(annotation, select=c("type"))

Var1        <- c("lightblue", "green")
names(Var1) <- c("tumor", "normal")
anno_colors <- list(type = Var1)

#my.cols <- greenred(100) #brewer.pal(9, "Blues")
pheatmap(scale(normalized_cds_sig), col = my.cols, annotation = annotation, cellwidth=3, cellheight=6,annotation_colors = anno_colors,
         show_colnames=F, legend=T, show_rownames=F, filename=paste(cancer, "/dz_sig_validation_lincs_reduced_v0.pdf", sep="")
)
#for some reason, in breast cancer, MIF has NA value after scale.
a = t(scale(t(normalized_cds_sig)))
a = a[!rownames(a) %in% c("MIF"), ]
pheatmap(a, col = my.cols, annotation = annotation, annotation_colors = anno_colors,
         show_colnames=F, legend=T, show_rownames=F, filename=paste(cancer, "/dz_sig_validation_lincs_reduced_v1.pdf", sep="")
)

