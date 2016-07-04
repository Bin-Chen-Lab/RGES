#compute correlation between RGES and drug IC50s after adding validation data of four compounds

library(cowplot)

###########
#functions##
###########
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
  return(sRGES * cor/max_cor) #
}
################

cancer = "LIHC"
landmark = 1
validated_drugs =  c("CGK-733", "strophanthidin", "FCCP",  "pyrvinium-pamoate")
strophanthidin = c(10, 0.9, 0.5, 1.8, 1.0)
FCCP = c(9.8, 1.9, 0.5, 1.6, 4.1)
CGK = c(3.4, 3.1, 2.1, 5.4, 2.1)
pyrvinium = c(0.07, 0.04, 0.04) 

validated_drugs_ic50 = c(median(CGK), median(strophanthidin), median(FCCP), median(pyrvinium)) * 1000

#build a reference model according to dose and time
output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
lincs_drug_prediction = read.csv(output_path)

lincs_drug_prediction_subset = subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
#pairs that share the same drug and cell id
lincs_drug_prediction_pairs = merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#x is the reference
lincs_drug_prediction_pairs = subset(lincs_drug_prediction_pairs, id.x != id.y & pert_time.x == 24 & pert_dose.x == 10) #, select = c("cmap_score.x", "cmap_score.y", "pert_dose.y", "pert_time.y"))

#difference of RGES to the reference 
lincs_drug_prediction_pairs$cmap_diff = lincs_drug_prediction_pairs$cmap_score.x - lincs_drug_prediction_pairs$cmap_score.y
lincs_drug_prediction_pairs$dose = round(log(lincs_drug_prediction_pairs$pert_dose.y, 2), 1)

#estimate difference
lincs_drug_prediction_pairs$dose_bin = ifelse(lincs_drug_prediction_pairs$pert_dose.y < 10, "low", "high")
tapply(lincs_drug_prediction_pairs$cmap_diff, lincs_drug_prediction_pairs$dose_bin, mean)
tapply(lincs_drug_prediction_pairs$cmap_diff, lincs_drug_prediction_pairs$pert_time.y, mean)
diff = tapply(lincs_drug_prediction_pairs$cmap_diff, paste(lincs_drug_prediction_pairs$dose_bin, lincs_drug_prediction_pairs$pert_time.y), mean)


#CMAP score output
cell_lines = read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))

#output_path <- paste(cancer, "/lincs_score_", landmark, ".csv", sep="")
output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
lincs_drug_prediction = read.csv(output_path)
lincs_drug_prediction$RGES = lincs_drug_prediction$cmap_score
lincs_drug_prediction = subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))

lincs_drug_activity = read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity = subset(lincs_drug_activity, standard_type == "IC50" &  cell_line %in% cell_lines$ChEMBL)
lincs_drug_activity = unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",		"cell_line")))
lincs_drug_activity = aggregate(standard_value ~ pert_iname , lincs_drug_activity, median)

lincs_drug_activity = rbind(lincs_drug_activity, data.frame(pert_iname= validated_drugs, standard_value=validated_drugs_ic50))
activity_RGES = merge(lincs_drug_prediction, lincs_drug_activity, by="pert_iname")

##weight cell lines
ccle_lincs = read.csv("raw/cell_line_lincs_ccle.csv")
cell_line_cancer = read.csv(paste(cancer, "/", "cell_line_", cancer, "_tacle.csv", sep=""))
cell_line_cancer$cor = cell_line_cancer$cor
cell_line_cancer = merge(cell_line_cancer, ccle_lincs, by.x="Cell.line.primary.name", by.y="ccle_cell_line_name")
cell_line_cancer = cell_line_cancer[order(cell_line_cancer$cor),]

results_subset = merge(activity_RGES, cell_line_cancer, by.x="cell_id", by.y="lincs_cell_id")

results_subset$RGES3 = sapply(1:nrow(results_subset), function(id){getsRGES3(results_subset$RGES[id], results_subset$cor[id], results_subset$pert_dose[id], results_subset$pert_time[id], diff, max(results_subset$cor))})

results_subset$sRGES = results_subset$RGES3

activity_RGES_summarized = aggregate(cbind(sRGES, standard_value) ~ pert_iname , results_subset, mean)

activity_RGES_summarized$p_text = sapply(1:nrow(activity_RGES_summarized), function(id){
  if (as.character(activity_RGES_summarized$pert_iname[id]) %in% validated_drugs){
    "red" #as.character(activity_RGES_summarized$pert_iname[id])
  }else{
    "black"
  }
})

pdf(paste( "fig/", cancer, "rges_ic50_normalized_validated_drugs.pdf", sep=""))
ggplot(activity_RGES_summarized, aes(sRGES, log(activity_RGES_summarized$standard_value, 10) )) +    theme_light() + theme_classic() + 
  theme(legend.position ="bottom", axis.text=element_text(size=22), axis.title=element_text(size=25)) +                                                                                               
  stat_smooth(method="lm", se=F, color="black")  + geom_point(aes(shape=p_text), size=3) + 
  scale_size(range = c(2, 5)) +
  xlab("sRGES") + guides(shape=FALSE, size=FALSE) +
  ylab("log10(IC50) nm") + coord_cartesian(xlim = c(-0.8, 0.8), ylim=c(-1, 8)) + 
  geom_vline(xintercept = -0.8) + geom_hline(yintercept = -1)

dev.off()

