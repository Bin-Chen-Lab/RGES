#compute correlation between RGES and drug IC50s; RGES is taken from the LINCS cloud
#first, need to use disease signatures to query LINCS and download the summary table and the detailed score table

library(cowplot)
library(pheatmap)

###########
#functions
##########
getsRGES1 <- function(RGES, cor, pert_dose, pert_time){
  sRGES <- RGES
  if (pert_time == 24){
    sRGES <- RGES + predict(lm_dose_24, data.frame(dose=round(log(pert_dose, 2), 1)))
  }
  if (pert_time == 6){
    sRGES <- RGES + predict(lm_dose_6, data.frame(dose=round(log(pert_dose, 2), 1)))
  }
  return (sRGES * cor )
}

getsRGES2 <- function(RGES, cor, pert_dose, pert_time){
  sRGES <- RGES 
  
  #older version
  if (pert_time < 24){
    sRGES <- sRGES - 0.1
  }
  
  if (pert_dose < 10){
    sRGES <- sRGES - 0.2
  }
  return(sRGES * cor)
}

getsRGES3 <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
  sRGES <- RGES
  pert_time <- ifelse(pert_time < 24, "short", "long")
  pert_dose <- ifelse(pert_dose < 10, "low", "high")
  if (pert_time == "short" & pert_dose == "low"){
    sRGES <- sRGES + diff[4]
  }
  if (pert_dose ==  "low" & pert_time == "long"){
    sRGES <- sRGES + diff[2]
  }
  if (pert_dose ==  "high" & pert_time == "short"){
    sRGES <- sRGES + diff[1]
  }
  return(sRGES * cor/max_cor) #* 
}

#####################
#MAIN
####################
cancer <- "BRCA"
type  <- "all" #all

srges_output_path <- paste(cancer, "/lincs_cancer_sRGES.csv", sep="")
summary_output_path <- paste(cancer, "/lincs_pred_summary_", type, ".txt", sep="")
raw_output_path <- paste(cancer, "/lincs_pred_all_", type, ".txt", sep="")

#preprocess data from the lincs cloud
lincs_drug_prediction <- read.csv(raw_output_path, sep="\t", stringsAsFactors = F)
lincs_drug_prediction <- subset(lincs_drug_prediction, pert_type == "trt_cp")
lincs_drug_prediction$pert_dose <- sapply(lincs_drug_prediction$pert_idose, function(x){
  items <- unlist(strsplit(x, " "))
  if (length(items) != 2){
    NA
  }else if (items[2] == "nM"){
    as.numeric(items[1])/1000
  }else{
    as.numeric(items[1])
  }})
lincs_drug_prediction$pert_dose_unit <- sapply(lincs_drug_prediction$pert_idose, function(x){unlist(strsplit(x, " "))[2]})
lincs_drug_prediction$pert_time <- as.numeric(sapply(lincs_drug_prediction$pert_itime, function(x){unlist(strsplit(x, " "))[1]}))
lincs_drug_prediction$pert_time_unit <- sapply(lincs_drug_prediction$pert_itime, function(x){unlist(strsplit(x, " "))[2]})
lincs_drug_prediction$cmap_score <- lincs_drug_prediction$connectivity_score
lincs_drug_prediction <- lincs_drug_prediction[!is.na(lincs_drug_prediction$pert_dose), ]
#

lincs_drug_prediction$id <- seq(1:nrow(lincs_drug_prediction))

lincs_drug_prediction_subset <- subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
#pairs that share the same drug and cell id
lincs_drug_prediction_pairs <- merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#x is the reference
lincs_drug_prediction_pairs <- subset(lincs_drug_prediction_pairs, id.x != id.y & pert_time.x == 24 & pert_dose.x == 10) #, select = c("cmap_score.x", "cmap_score.y", "pert_dose.y", "pert_time.y"))

#difference of RGES to the reference 
lincs_drug_prediction_pairs$cmap_diff <- lincs_drug_prediction_pairs$cmap_score.x - lincs_drug_prediction_pairs$cmap_score.y
lincs_drug_prediction_pairs$dose <- round(log(lincs_drug_prediction_pairs$pert_dose.y, 2), 1)

#fix time
lincs_drug_prediction_pairs_subset <- subset(lincs_drug_prediction_pairs, pert_time.y == 24 )
dose_cmap_diff_24 <- tapply(lincs_drug_prediction_pairs_subset$cmap_diff, lincs_drug_prediction_pairs_subset$dose, mean)
dose_cmap_diff_24 <- data.frame(dose = as.numeric(names(dose_cmap_diff_24)), cmap_diff= dose_cmap_diff_24)
plot(dose_cmap_diff_24$dose, dose_cmap_diff_24$cmap_diff)
lm_dose_24 <- lm(cmap_diff ~ dose, data = dose_cmap_diff_24)
summary(lm_dose_24)

lincs_drug_prediction_pairs_subset <- subset(lincs_drug_prediction_pairs, pert_time.y == 6)
dose_cmap_diff_6 <- tapply(lincs_drug_prediction_pairs_subset$cmap_diff, lincs_drug_prediction_pairs_subset$dose, mean)
dose_cmap_diff_6 <- data.frame(dose = as.numeric(names(dose_cmap_diff_6)), cmap_diff= dose_cmap_diff_6)
lm_dose_6 <- lm(cmap_diff ~ dose, data = dose_cmap_diff_6)
plot(dose_cmap_diff_6$dose, dose_cmap_diff_6$cmap_diff)
summary(lm_dose_6)

#estimate difference
lincs_drug_prediction_pairs$dose_bin <- ifelse(lincs_drug_prediction_pairs$pert_dose.y < 10, "low", "high")
tapply(lincs_drug_prediction_pairs$cmap_diff, lincs_drug_prediction_pairs$dose_bin, mean)
tapply(lincs_drug_prediction_pairs$cmap_diff, lincs_drug_prediction_pairs$pert_time.y, mean)
diff <- tapply(lincs_drug_prediction_pairs$cmap_diff, paste(lincs_drug_prediction_pairs$dose_bin, lincs_drug_prediction_pairs$pert_time.y), mean)


#CMAP score output
cell_lines <- read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))

lincs_drug_prediction$RGES <- lincs_drug_prediction$connectivity_score
lincs_drug_prediction <- subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))

lincs_drug_activity <- read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity <- subset(lincs_drug_activity, standard_type == "IC50" &  cell_line %in% cell_lines$ChEMBL)
lincs_drug_activity <- unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",		"cell_line")))
lincs_drug_activity <- aggregate(standard_value ~ pert_iname , lincs_drug_activity, median)

activity_RGES <- merge(lincs_drug_prediction, lincs_drug_activity, by="pert_iname")

###same lineage
activity_RGES_subset <- subset(activity_RGES, (cell_id %in% cell_lines$LINCS)) 
activity_RGES_summarized <- aggregate(cbind(RGES, standard_value) ~ pert_iname, activity_RGES_subset,  median)
cor_test <- cor.test(activity_RGES_summarized$RGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test

activity_RGES_summarized <- aggregate(cbind(RGES, standard_value) ~ pert_iname, activity_RGES,  median)
t.test(activity_RGES_summarized$RGES[activity_RGES_summarized$standard_value < 10000], activity_RGES_summarized$RGES[activity_RGES_summarized$standard_value > 10000])

cor_test <- cor.test(activity_RGES_summarized$RGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test

lm_cmap_ic50 <- lm(RGES ~ log(standard_value, 10), activity_RGES_summarized)
summary(lm_cmap_ic50)
summary(lm_cmap_ic50)$r.squared

pdf(paste( "fig/", cancer, "rges_ic50_min_lincs_cloud_", type, ".pdf", sep=""))
ggplot(activity_RGES_summarized, aes(RGES, log(activity_RGES_summarized$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label <- paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=2), sep=""), x = 0, y = 7.7, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 7.3, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("RGES") + guides(shape=FALSE, size=FALSE) +
  ylab("log10(IC50) nm") + coord_cartesian(xlim = c(-1, 1), ylim=c(-1, 8)) 
dev.off()


##weight cell lines
ccle_lincs <- read.csv("raw/cell_line_lincs_ccle.csv")
cell_line_cancer <- read.csv(paste( cancer,  "/", "cell_line_", cancer, "_tacle.csv", sep=""))
cell_line_cancer <- merge(cell_line_cancer, ccle_lincs, by.x="Cell.line.primary.name", by.y="ccle_cell_line_name")
cell_line_cancer <- cell_line_cancer[order(cell_line_cancer$cor),]

results_subset <- merge(activity_RGES, cell_line_cancer, by.x="cell_id", by.y="lincs_cell_id")

#should keep the cell lines from the same lineage
#results_subset <- subset(results_subset, !(cell_id %in% cell_lines$LINCS))
results_subset$RGES1 <- sapply(1:nrow(results_subset), function(id){getsRGES1(results_subset$RGES[id], results_subset$cor[id], results_subset$pert_dose[id], results_subset$pert_time[id])})
results_subset$RGES2 <- sapply(1:nrow(results_subset), function(id){getsRGES2(results_subset$RGES[id], results_subset$cor[id], results_subset$pert_dose[id], results_subset$pert_time[id])})
results_subset$RGES3 <- sapply(1:nrow(results_subset), function(id){getsRGES3(results_subset$RGES[id], results_subset$cor[id], results_subset$pert_dose[id], results_subset$pert_time[id], diff, max(results_subset$cor))})

results_subset$sRGES <- results_subset$RGES3

activity_RGES_summarized <- aggregate(cbind(sRGES, standard_value) ~ pert_iname , results_subset, mean)

lincs_drug_prediction_rges <- activity_RGES_summarized

cor_test <- cor.test(activity_RGES_summarized$sRGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test

lm_cmap_ic50 <- lm( log(standard_value, 10) ~ sRGES, activity_RGES_summarized)
summary(lm_cmap_ic50)
summary(lm_cmap_ic50)$r.squared^0.5


activity_RGES_summarized$activity = "Active (<=4)"
activity_RGES_summarized$activity[activity_RGES_summarized$standard_value>10000] = "Inactive (>4)"

pdf(paste( "fig/", cancer, "rges_ic50_normalized_lincs_cloud_", type, ".pdf", sep=""))
ggplot(activity_RGES_summarized, aes(sRGES, log(activity_RGES_summarized$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18)) +                                                                                               
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label =  cancer , x = 0, y = 7.9, size = 7) + 
  annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=3, scientific=T), sep=""), x = 0, y = 7.5, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 7.1, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("sRGES") + guides(shape=FALSE, size=FALSE) +
  ylab("log10(IC50) nm") + coord_cartesian(xlim = c(-0.7, 0.7), ylim=c(-1, 8))
dev.off()

############################
##compare with our sRGES from the summarization method from the lincs cloud
lincs_drug_prediction_summary <- read.csv(summary_output_path, sep="\t", stringsAsFactors = F)
results_subset <- merge(lincs_drug_prediction_summary, lincs_drug_activity, by = "pert_iname")

results_subset$sRGES <- results_subset$mean_rankpt_4 #mean_rankpt_4 was recommended in the LINCS cloud
activity_RGES_summarized <- aggregate(cbind(sRGES, standard_value) ~ pert_iname , results_subset, mean)
cor_test <- cor.test(activity_RGES_summarized$sRGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test

lm_cmap_ic50 <- lm( log(standard_value, 10) ~ sRGES, activity_RGES_summarized)
summary(lm_cmap_ic50)
summary(lm_cmap_ic50)$r.squared^0.5

pdf(paste( "fig/", cancer, "rges_ic50_sRGES_vs_lincs_summarization_", type, ".pdf", sep=""))
ggplot(activity_RGES_summarized, aes(sRGES, log(activity_RGES_summarized$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18)) +                                                                                               
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label =  cancer , x = 0, y = 7.9, size = 7) + 
  annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=3, scientific=T), sep=""), x = 0, y = 7.5, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 7.1, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("sRGES") + guides(shape=FALSE, size=FALSE) +
  ylab("log10(IC50) nm") + coord_cartesian(xlim = c(-100, 100), ylim=c(-1, 8))
dev.off()

#  
###compare sRGES with the summarization method from the LINCS cloud
pred_all <- merge(lincs_drug_prediction_rges, lincs_drug_prediction_summary, by="pert_iname")
pred_all$IC50 <- pred_all$standard_value
pred_all$summarized_CMap_score <- pred_all$sRGES
pred_all <- pred_all[, c("pert_iname",  "IC50", "summarized_CMap_score", "mean_rankpt_2", "mean_rankpt_4", "mean_rankpt_6", "rankpt_1", "rankpt_2", "rankpt_3",
                        "rankpt_4", "rankpt_5", "rankpt_6", "rankpt_7", "rankpt_9")]
pred_all_merged <- aggregate(. ~ pert_iname, pred_all, mean)
#View(pred_all_merged)
pred_all_merged <- pred_all_merged[, -1]
pred_all_merged_cor <- cor(pred_all_merged, method="spearman")
pred_all_merged_cor["IC50",]

my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = 50)
pheatmap(pred_all_merged_cor, color = my_palette, cluster_rows=F, cluster_cols=F, 
         cellwidth=25,  cellheight=25, display_numbers=T, file=paste("fig/", cancer, "_", type, "_performance_with_srges.pdf", sep=""))

write.csv(pred_all_merged_cor, paste(cancer, "/compare_sRGES_lincs_cloud_", type, ".csv", sep=""))

###compare with the results from our method using our LINCS data library
lincs_drug_prediction_rges <- read.csv(srges_output_path,  stringsAsFactors = F)
pred_all <- merge(lincs_drug_prediction_rges, lincs_drug_prediction_summary, by="pert_iname")
pred_all <- merge(pred_all, lincs_drug_activity, by="pert_iname")
pred_all$sRGES <- pred_all$sRGES
pred_all$IC50 <- pred_all$standard_value
pred_all <- pred_all[, c("pert_iname",  "IC50", "sRGES", "mean_rankpt_2", "mean_rankpt_4", "mean_rankpt_6", "rankpt_1", "rankpt_2", "rankpt_3",
                        "rankpt_4", "rankpt_5", "rankpt_6", "rankpt_7", "rankpt_9")]
pred_all_merged <- aggregate(. ~ pert_iname, pred_all, mean)
#View(pred_all_merged)
pred_all_merged <- pred_all_merged[, -1]
pred_all_merged_cor <- cor(pred_all_merged, method="spearman")
pred_all_merged_cor["IC50",]

my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = 50)
pheatmap(pred_all_merged_cor, color = my_palette, cluster_rows=F, cluster_cols=F, 
         cellwidth=25,  cellheight=25, display_numbers=T, file=paste("fig/", cancer, "_", type, "_performance.pdf", sep=""))



