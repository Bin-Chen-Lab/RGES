#Compute RGES in CMap and evaluate if the correlation is retained
#only BRCA has signficant correlation; LIHC and COAD do not have related cell lines.
#correlate drug sensitivity with the potency to reserse disease gene expression (RGES)
#to reduce the computation time, only the drugs with sensitivty data were examined

cancer <- "BRCA"
landmark <- 0
#CMAP score output
output_path <- paste(cancer, "/lincs_score_", "geo", ".csv", sep="")

##############
cmap_score_new <- function(sig_up, sig_down, drug_signature) {

  num_genes <- nrow(drug_signature)
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

get.gene.list <- function(con){
  lincs_gene <- dbReadTable(con, "probe_id_info")
  
  return(lincs_gene$gene_id)
}


##############
load("raw/lincs/lincs_signatures_cmpd_landmark.RData")
load('raw/cmap//cmap_signatures.RData')
cmap_sig_info <- read.csv("raw/cmap/cmap_drug_experiments_new.csv", stringsAsFactors=F)
cmap_lincs_valid_instances <- read.csv("raw/cmap/cmap_valid_instances.csv")

#cmap_sig_info <- subset(cmap_sig_info, cell_line %in% c("MCF7") & concentration < 5*1E-5 & concentration> 5*1E-6 & id %in% cmap_lincs_valid_instances$id[cmap_lincs_valid_instances$valid %in% c(1)])
cmap_sig_info <- subset(cmap_sig_info, cell_line %in% c("MCF7") & id %in% cmap_lincs_valid_instances$id[cmap_lincs_valid_instances$valid %in% c(1)])

lincs_drug_activity <- read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity <- unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",		"cell_line")))

landmark <- 0
if (landmark ==1){
  cmap_signatures <- cmap_signatures[cmap_signatures$V1 %in% rownames(lincs_signatures), ]
}
gene.list <- cmap_signatures$V1

sig.ids <- cmap_sig_info$id[tolower(cmap_sig_info$name) %in% tolower(lincs_drug_activity$pert_iname)]

#load  genes
fromGSE62944 <- F
if (fromGSE62944){
  dz_signature <- read.table(dz_sig_path,header=T,sep="\t")
  dz_signature <- subset(dz_signature, GeneID %in% gene.list)
}else{
  load(paste(cancer, '/', cancer, 'tcga_deseq_final.RData', sep='')) #tcga_deseq LIHC_sig_deseq_QC
  res$GeneID <- sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[2]})
  res$symbol <- sapply((res$id), function(id){unlist(strsplit(id, "\\|"))[1]})
  dz_signature <- subset(res, !is.na(padj) & !is.na(id) & id !='?' & padj < 1E-3 & abs(log2FoldChange) > 1.5 & abs(log2FoldChange) != Inf )
  dz_signature <- subset(dz_signature, GeneID %in% gene.list ) #& GeneID %in% rownames(lincs_signatures))
  dz_signature$up_down <- "up"
  dz_signature$up_down[dz_signature$log2FoldChange<0] <- "down"
}

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
  sig <-  cmap_signatures[, paste("V", as.character(exp_id+1), sep="")]
  cmap_exp_signature <- data.frame(gene.list,  rank(-1 * sig, ties.method="random"))    
  colnames(cmap_exp_signature) <- c("ids","rank")
  dz_cmap_scores <- c(dz_cmap_scores, cmap_score_new(dz_genes_up,dz_genes_down,cmap_exp_signature))
  
 # compute spearman correlation and pearson correlation
  dz_sig_drug_sig <- merge(dz_signature, data.frame(GeneID=gene.list, expr= sig), by="GeneID")
  dz_spearmans <- c(dz_spearmans, cor(dz_sig_drug_sig$log2FoldChange, dz_sig_drug_sig$expr, method="spearman"))
  IQRs <- c(IQRs, IQR(sig))
  
}

results <- data.frame(id = sig.ids, RGES = dz_cmap_scores, spearman = dz_spearmans, iqr = IQRs)

results <- merge(results, cmap_sig_info, by = "id")
results <- results[order(results$RGES),]
write.csv(results, output_path)

##########
#correlated to IC50
lincs_drug_activity <- aggregate(standard_value ~ pert_iname, lincs_drug_activity,median)

drug_activity_rges <- merge(results, lincs_drug_activity, by.x="name", by.y="pert_iname")

drug_activity_rges <- aggregate(cbind(RGES, standard_value) ~ name, drug_activity_rges, median)

plot(drug_activity_rges$RGES, log(drug_activity_rges$standard_value, 10))
cor_test <- cor.test(drug_activity_rges$RGES, log(drug_activity_rges$standard_value, 10))

drug_activity_rges <- drug_activity_rges[order(drug_activity_rges$RGES),]

lm_cmap_ic50 <- lm(RGES ~ log(standard_value, 10), drug_activity_rges)


pdf(paste( "fig/", cancer, "rges_ic50_cmap_data_", landmark, ".pdf", sep=""))
ggplot(drug_activity_rges, aes(RGES, log(drug_activity_rges$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label = paste(cancer, ",", "MCF7", sep=""), 
           x = 0, y = 8.1, size = 6, colour = "black") +
  annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=2), sep=""), 
           x = 0, y = 7.7, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 7.3, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("RGES") + guides(shape=FALSE, size=FALSE) +
  ylab("log10(IC50) nm") + coord_cartesian(xlim = c(-0.5, 0.5), ylim=c(-1, 8)) 
dev.off()

