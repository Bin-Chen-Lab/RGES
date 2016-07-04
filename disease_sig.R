#build up disease signatures

cancer = "ER"

load("raw/lincs/lincs_signatures_cmpd_landmark.RData")
gene.list = rownames(lincs_signatures)

load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC
res$GeneID = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
res$symbol = sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
dz_signature = subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf )

dz_signature$landmark = sign(match(dz_signature$GeneID,  gene.list ))

a1 = dz_signature[!is.na(dz_signature$landmark), ]
write.csv(dz_signature, paste(cancer, '/', cancer, "_sig.csv", sep=""))

