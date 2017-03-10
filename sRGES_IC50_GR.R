#compute correlation between RGES and drug IC50s

library(cowplot)
library("ROCR")
library(pROC)
library(ggplot2)
##################
#functions
##################
getsRGES1 = function(RGES, cor, pert_dose, pert_time){
  sRGES = RGES
  if (pert_time == 24){
    sRGES = RGES + predict(lm_dose_24, data.frame(dose=round(log(pert_dose, 10), 1)))
  }
  if (pert_time == 6){
    sRGES = RGES + predict(lm_dose_6, data.frame(dose=round(log(pert_dose, 10), 1)))
  }
  return (sRGES * cor )
}

getsRGES2 = function(RGES, cor, pert_dose, pert_time){
  sRGES = RGES 
  
  #older version
  if (pert_time < 24){
    sRGES = sRGES - 0.1
  }
  
  if (pert_dose < 10){
    sRGES = sRGES - 0.2
  }
  return(sRGES * cor)
}

getsRGES3 = function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
  sRGES = RGES
  pert_time = ifelse(pert_time < 24, "short", "long")
  pert_dose = ifelse(pert_dose < 10, "low", "high")
  if (pert_time == "short" & pert_dose == "low"){
    sRGES = sRGES + diff[4]
  }
  if (pert_dose ==  "low" & pert_time == "long"){
    sRGES = sRGES + diff[2]
  }
  if (pert_dose ==  "high" & pert_time == "short"){
    sRGES = sRGES + diff[1]
  }
  return(sRGES ) #* cor/max_cor
}


#################
#MAIN
################

cancer = "BRCA"
landmark = 1 #landmark: 1. all genes: 0

#build a reference model according to dose and time
output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
lincs_drug_prediction = read.csv(output_path)

#should use pert_dose > 0.01?
lincs_drug_prediction_subset = subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
#pairs that share the same drug and cell id
lincs_drug_prediction_pairs = merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#x is the reference
lincs_drug_prediction_pairs = subset(lincs_drug_prediction_pairs, id.x != id.y & pert_time.x == 24 & pert_dose.x == 10) #, select = c("cmap_score.x", "cmap_score.y", "pert_dose.y", "pert_time.y"))

#difference of RGES to the reference 
lincs_drug_prediction_pairs$cmap_diff = lincs_drug_prediction_pairs$cmap_score.x - lincs_drug_prediction_pairs$cmap_score.y
lincs_drug_prediction_pairs$dose = round(log(lincs_drug_prediction_pairs$pert_dose.y, 10), 1)

#fix time
lincs_drug_prediction_pairs_subset = subset(lincs_drug_prediction_pairs, pert_time.y == 24 )
dose_cmap_diff_24 = tapply(lincs_drug_prediction_pairs_subset$cmap_diff, lincs_drug_prediction_pairs_subset$dose, mean)
dose_cmap_diff_24 = data.frame(dose = as.numeric(names(dose_cmap_diff_24)), cmap_diff= dose_cmap_diff_24)
plot(dose_cmap_diff_24$dose, dose_cmap_diff_24$cmap_diff)
lm_dose_24 = lm(cmap_diff ~ dose, data = dose_cmap_diff_24)
summary(lm_dose_24)

lincs_drug_prediction_pairs_subset = subset(lincs_drug_prediction_pairs, pert_time.y == 6)
dose_cmap_diff_6 = tapply(lincs_drug_prediction_pairs_subset$cmap_diff, lincs_drug_prediction_pairs_subset$dose, mean)
dose_cmap_diff_6 = data.frame(dose = as.numeric(names(dose_cmap_diff_6)), cmap_diff= dose_cmap_diff_6)
lm_dose_6 = lm(cmap_diff ~ dose, data = dose_cmap_diff_6)
plot(dose_cmap_diff_6$dose, dose_cmap_diff_6$cmap_diff)
summary(lm_dose_6)

#estimate difference
lincs_drug_prediction_pairs$dose_bin = ifelse(lincs_drug_prediction_pairs$pert_dose.y < 10, "low", "high")
tapply(lincs_drug_prediction_pairs$cmap_diff, lincs_drug_prediction_pairs$dose_bin, mean)
tapply(lincs_drug_prediction_pairs$cmap_diff, lincs_drug_prediction_pairs$pert_time.y, mean)
diff = tapply(lincs_drug_prediction_pairs$cmap_diff, paste(lincs_drug_prediction_pairs$dose_bin, lincs_drug_prediction_pairs$pert_time.y), mean)


#CMAP score output
cell_lines = read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))
if ("subtype" %in% colnames(cell_lines)){ cell_lines = cell_lines[!is.na(cell_lines$subtype),]}

#output_path <- paste(cancer, "/lincs_score_", landmark, ".csv", sep="")
output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
lincs_drug_prediction = read.csv(output_path)
lincs_drug_prediction$RGES = lincs_drug_prediction$cmap_score
lincs_drug_prediction = subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))

lincs_drug_activity = read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity = subset(lincs_drug_activity, standard_type == "IC50" &  cell_line %in% cell_lines$ChEMBL) #
lincs_drug_activity = unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",		"cell_line")))
lincs_drug_activity = aggregate(standard_value ~ pert_iname , lincs_drug_activity, median)

#data_4_Heiser_2012.csv
#data_1_LINCS_Pilot_Phase_Joint_Project.csv
lincs_drug_activity_gr = read.csv(paste(cancer, "/data_1_LINCS_Pilot_Phase_Joint_Project.csv", sep=""), stringsAsFactors = F)
lincs_drug_activity_gr$pert_iname = tolower(lincs_drug_activity_gr$smallMolecule)
lincs_drug_activity_gr$standard_value = as.numeric(lincs_drug_activity_gr$GRmax)
#lincs_drug_activity_gr$cellLine %in% c("MCF7") & 
lincs_drug_activity_gr = lincs_drug_activity_gr[ !is.na(lincs_drug_activity_gr$standard_value) & lincs_drug_activity_gr$standard_value != "Inf",]
lincs_drug_activity_gr = aggregate(standard_value ~ pert_iname, lincs_drug_activity_gr, median)
lincs_drug_activity_gr$standard_value = 10 ^ lincs_drug_activity_gr$standard_value

ic50_gr = merge(lincs_drug_activity, lincs_drug_activity_gr, by = "pert_iname")
cor.test(ic50_gr$standard_value.x, ic50_gr$standard_value.y, method="spearman")
plot(log(ic50_gr$standard_value.x), ic50_gr$standard_value.y)

activity_RGES = merge(lincs_drug_prediction, lincs_drug_activity_gr, by="pert_iname")

###same lineage
activity_RGES_subset = subset(activity_RGES, (cell_id %in% cell_lines$LINCS)) 
activity_RGES_summarized = aggregate(cbind(RGES, standard_value) ~ pert_iname, activity_RGES_subset,  median)
cor_test = cor.test(activity_RGES_summarized$RGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test

#igoring normal_cells may increase the correlation, tweaking dose and time does not. 
#ignoring neither blood cell or tnbc cells may increase the correlation
activity_RGES_summarized = aggregate(cbind(RGES, standard_value) ~ pert_iname, activity_RGES,  mean)
t.test(activity_RGES_summarized$RGES[activity_RGES_summarized$standard_value < 10000], activity_RGES_summarized$RGES[activity_RGES_summarized$standard_value > 10000])


lm_cmap_ic50 = lm(RGES ~ log(standard_value, 10), activity_RGES_summarized)
summary(lm_cmap_ic50)
summary(lm_cmap_ic50)$r.squared

activity_RGES$activity = 1
activity_RGES$activity[activity_RGES$standard_value>10000] = 2

##weight cell lines
ccle_lincs = read.csv("raw/cell_line_lincs_ccle.csv")
cell_line_cancer = read.csv(paste(cancer, "/", "cell_line_", cancer, "_tacle.csv", sep=""))
cell_line_cancer = merge(cell_line_cancer, ccle_lincs, by.x="Cell.line.primary.name", by.y="ccle_cell_line_name")
cell_line_cancer = cell_line_cancer[order(cell_line_cancer$cor),]

results_subset = merge(activity_RGES, cell_line_cancer, by.x="cell_id", by.y="lincs_cell_id")

#should keep the cell lines from the same lineage
#results_subset = subset(results_subset, !(cell_id %in% cell_lines$LINCS))

results_subset$RGES1 = sapply(1:nrow(results_subset), function(id){getsRGES1(results_subset$RGES[id], results_subset$cor[id], results_subset$pert_dose[id], results_subset$pert_time[id])})
results_subset$RGES2 = sapply(1:nrow(results_subset), function(id){getsRGES2(results_subset$RGES[id], results_subset$cor[id], results_subset$pert_dose[id], results_subset$pert_time[id])})
results_subset$RGES3 = sapply(1:nrow(results_subset), function(id){getsRGES3(results_subset$RGES[id], results_subset$cor[id], results_subset$pert_dose[id], results_subset$pert_time[id], diff, max(results_subset$cor))})

#using RGES3 by default
results_subset$sRGES = results_subset$RGES3

activity_RGES_summarized = aggregate(cbind(sRGES, standard_value) ~ pert_iname , results_subset, mean)

cor_test = cor.test(activity_RGES_summarized$sRGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test
lm_cmap_ic50 = lm( log(standard_value, 10) ~ sRGES, activity_RGES_summarized)
summary(lm_cmap_ic50)
summary(lm_cmap_ic50)$r.squared^0.5


activity_RGES_summarized$activity = "Active (<=4)"
activity_RGES_summarized$activity[activity_RGES_summarized$standard_value>10000] = "Inactive (>4)"

pdf(paste( "fig/", cancer, "rges_gr_normalized.pdf", sep=""))
ggplot(activity_RGES_summarized, aes(sRGES, log(activity_RGES_summarized$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18)) +                                                                                               
 stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
      annotate("text", label =  cancer , x = 0, y = 2, size = 7) + 
        annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=3, scientific=T), sep=""), x = 0, y = 1.7, size = 6, colour = "black") +
        annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 1.6, size = 6, colour = "black") +
        scale_size(range = c(2, 5)) +
        xlab("sRGES") + guides(shape=FALSE, size=FALSE) +
        ylab("GR max") + coord_cartesian(xlim = c(-0.7, 0.7), ylim=c(-1, 2))
dev.off()

