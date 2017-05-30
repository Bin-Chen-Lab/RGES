#an example of running RGES and summarizing RGES across multiple profiles.
setwd("~/Documents/stanford/tumor_cell_line/RGES_manuscript/release/data")

library("plyr")
library("ggplot2")

load("raw/lincs/lincs_signatures_cmpd_landmark.RData")

#
code_dir <- "../git_code/RGES/"
source(paste(code_dir, "core_functions.R",sep=""))

output_path <- paste("example/all_lincs_score.csv", sep="")
sRGES_output_path <- paste("example/sRGES.csv", sep="")

landmark <- 1
lincs_sig_info <- read.csv("raw/lincs/lincs_sig_info.csv")
lincs_sig_info <- subset(lincs_sig_info, id %in% colnames(lincs_signatures))
#remove duplicate instances
lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]

sig.ids <- lincs_sig_info$id

##############
#read disease signatures
dz_signature <- read.csv("example/LIHC_sig_reduced.csv")
dz_genes_up <- subset(dz_signature,up_down=="up",select="GeneID")
dz_genes_down <- subset(dz_signature,up_down=="down",select="GeneID")
###############

##############
#compute RGES
#only support landmark genes
gene.list <- rownames(lincs_signatures)

#compute RGES
#only choose the top 100 genes
max_gene_size <- 150
if (nrow(dz_genes_up) > max_gene_size){
  dz_genes_up <- data.frame(GeneID= dz_genes_up[1:max_gene_size,])
}
if (nrow(dz_genes_down) > max_gene_size){
  dz_genes_down <- data.frame(GeneID=dz_genes_down[1:max_gene_size,])
}

dz_cmap_scores <- NULL
count <- 0
for (exp_id in sig.ids) {
  count <- count + 1
  print(count)
  if (landmark ==1){
    cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
  }else{
    cmap_exp_signature <- cbind(gene.list,  get.sigs(exp_id))    
  }
  colnames(cmap_exp_signature) <- c("ids","rank")
  dz_cmap_scores <- c(dz_cmap_scores, cmap_score_new(dz_genes_up,dz_genes_down,cmap_exp_signature))
}

#random scores
N_PERMUTATIONS <- 10000 #default 100000
random_sig_ids <- sample(1:ncol(lincs_signatures),N_PERMUTATIONS,replace=T)
count <- 0
random_cmap_scores <- NULL
for (expr_id in random_sig_ids){
  count <- count + 1
  print(count)
  if (landmark ==1){
    cmap_exp_signature <- data.frame(gene.list,  rank(-1 * lincs_signatures[, as.character(exp_id)], ties.method="random"))    
  }else{
    cmap_exp_signature <- cbind(gene.list,  get.sigs(exp_id))    
  }
  colnames(cmap_exp_signature) <- c("ids","rank")
  
  random_input_signature_genes <- sample(gene.list, (nrow(dz_genes_up)+nrow(dz_genes_down)))
  rand_dz_gene_up <- data.frame(GeneID=random_input_signature_genes[1:nrow(dz_genes_up)])
  rand_dz_gene_down <- data.frame(GeneID=random_input_signature_genes[(nrow(dz_genes_up)+1):length(random_input_signature_genes)])
  random_cmap_scores <- c(random_cmap_scores, cmap_score_new(rand_dz_gene_up,rand_dz_gene_down,cmap_exp_signature))
}

p <- sapply(dz_cmap_scores, function(score){
  sum(random_cmap_scores < score)/length(random_cmap_scores)
})

padj <- p.adjust(p, "fdr")
results <- data.frame(id = sig.ids, cmap_score = dz_cmap_scores, p, padj)

results <- merge(results, lincs_sig_info, by = "id")
results <- results[order(results$cmap_score),]
write.csv(results, output_path)

####################
#summarize RGES
lincs_drug_prediction <- read.csv(output_path)

#should use pert_dose > 0.01
lincs_drug_prediction_subset <- subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
#pairs that share the same drug and cell id
lincs_drug_prediction_pairs <- merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#x is the reference
lincs_drug_prediction_pairs <- subset(lincs_drug_prediction_pairs, id.x != id.y & pert_time.x == 24 & pert_dose.x == 10) #, select <- c("cmap_score.x", "cmap_score.y", "pert_dose.y", "pert_time.y"))

#difference of RGES to the reference 
lincs_drug_prediction_pairs$cmap_diff <- lincs_drug_prediction_pairs$cmap_score.x - lincs_drug_prediction_pairs$cmap_score.y
lincs_drug_prediction_pairs$dose <- round(log(lincs_drug_prediction_pairs$pert_dose.y, 2), 1)

#estimate difference
lincs_drug_prediction_pairs$dose_bin <- ifelse(lincs_drug_prediction_pairs$pert_dose.y < 10, "low", "high")
diff <- tapply(lincs_drug_prediction_pairs$cmap_diff, paste(lincs_drug_prediction_pairs$dose_bin, lincs_drug_prediction_pairs$pert_time.y), mean)

#ignore weighting cell lines
lincs_cell_line_weight <- read.csv("example/lincs_cell_line_weight.csv")
pred <- merge(lincs_drug_prediction, lincs_cell_line_weight, by.x="cell_id", by.y="lincs_cell_id")
pred$RGES <- sapply(1:nrow(pred), function(id){getsRGES(pred$cmap_score[id], pred$cor[id], pred$pert_dose[id], pred$pert_time[id], diff, max(pred$cor))})

cmpd_freq <- table(pred$pert_iname)
pred <- subset(pred, pert_iname %in% names(cmpd_freq[cmpd_freq > 0]))

pred_merged <- ddply(pred,  .(pert_iname),  summarise,
                     mean = mean(RGES),
                     n = length(RGES),
                     median = median(RGES),
                     sd = sd(RGES))
pred_merged$sRGES <- pred_merged$mean
pred_merged <- pred_merged[order(pred_merged$sRGES), ]

write.csv(pred_merged, sRGES_output_path)

####
#sRGES and drug efficacy
drug_efficacy <- read.csv("example/drug_efficacy.csv")

RGES_efficacy <- merge(pred_merged, drug_efficacy, by = "pert_iname")

cor_test <- cor.test(RGES_efficacy$sRGES, log(RGES_efficacy$standard_value, 10), method="spearman", exact=FALSE) #or kendall
lm_cmap_ic50 <- lm( log(standard_value, 10) ~ sRGES, RGES_efficacy)

ggplot(RGES_efficacy, aes(sRGES, log(RGES_efficacy$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18)) +                                                                                               
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=3, scientific=T), sep=""), x = 0, y = 7.5, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 7.1, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("sRGES") + guides(shape=FALSE, size=FALSE) +
  ylab("log10(IC50) nm") + coord_cartesian(xlim = c(-0.7, 0.7), ylim=c(-1, 8))

