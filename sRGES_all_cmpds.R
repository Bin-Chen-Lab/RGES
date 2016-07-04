#using sRGES to summarize all compounds

library("plyr")
##################
####
getsRGES = function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
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

##############
##############

cancer = "ER"

#build a reference model according to dose and time
output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
lincs_drug_prediction = read.csv(output_path)

#should use pert_dose > 0.01
lincs_drug_prediction_subset = subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
#pairs that share the same drug and cell id
lincs_drug_prediction_pairs = merge(lincs_drug_prediction_subset, lincs_drug_prediction_subset, by=c("pert_iname", "cell_id")) 
#x is the reference
lincs_drug_prediction_pairs = subset(lincs_drug_prediction_pairs, id.x != id.y & pert_time.x == 24 & pert_dose.x == 10) #, select = c("cmap_score.x", "cmap_score.y", "pert_dose.y", "pert_time.y"))

#difference of RGES to the reference 
lincs_drug_prediction_pairs$cmap_diff = lincs_drug_prediction_pairs$cmap_score.x - lincs_drug_prediction_pairs$cmap_score.y
lincs_drug_prediction_pairs$dose = round(log(lincs_drug_prediction_pairs$pert_dose.y, 2), 1)

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


cell_lines = read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))

cell_line_cancer = read.csv(paste(cancer, "/", "cell_line_", cancer, "_tacle.csv", sep=""))
cell_line_cancer = merge(cell_line_cancer, ccle_lincs, by.x="Cell.line.primary.name", by.y="ccle_cell_line_name")
cell_line_cancer = cell_line_cancer[order(cell_line_cancer$cor),]

pred = merge(lincs_drug_prediction, cell_line_cancer, by.x="cell_id", by.y="lincs_cell_id")

pred$RGES = sapply(1:nrow(pred), function(id){getsRGES(pred$cmap_score[id], pred$cor[id], pred$pert_dose[id], pred$pert_time[id], diff, max(pred$cor))})

cmpd_freq = table(pred$pert_iname)
pred = subset(pred, pert_iname %in% names(cmpd_freq[cmpd_freq>0]))

pred_merged = ddply(pred,  .(pert_iname),  summarise,
                        mean = mean(RGES),
                        n = length(RGES),
                        median = median(RGES),
                        sd = sd(RGES))
pred_merged$sRGES = pred_merged$mean
pred_merged = pred_merged[order(pred_merged$sRGES), ]

write.csv(pred_merged,paste( cancer, "/lincs_cancer_sRGES.csv", sep=""))


