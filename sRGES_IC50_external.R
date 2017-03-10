#compute correlation between RGES and drug IC50s

library(cowplot)
library("ROCR")
library(pROC)

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

cancer = "ER"
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

output_path <- paste(cancer, "/lincs_score_", landmark, ".csv", sep="")
lincs_drug_prediction = read.csv(output_path)
lincs_drug_prediction$RGES = lincs_drug_prediction$RGES
lincs_drug_prediction = subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))

lincs_drug_activity = read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity = subset(lincs_drug_activity, standard_type == "IC50" &  cell_line %in% cell_lines$ChEMBL)
lincs_drug_activity = unique(subset(lincs_drug_activity, select=c("pert_iname", "doc_id", "standard_value", "standard_type", "description",    "organism",		"cell_line")))
lincs_drug_activity = aggregate(standard_value ~ pert_iname , lincs_drug_activity, median)

activity_RGES = merge(lincs_drug_prediction, lincs_drug_activity, by="pert_iname")

###same lineage
activity_RGES_subset = subset(activity_RGES, (cell_id %in% cell_lines$LINCS)) 
activity_RGES_summarized = aggregate(cbind(RGES, spearman, pearson, standard_value) ~ pert_iname, activity_RGES_subset,  median)
cor_test = cor.test(activity_RGES_summarized$RGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test

#igoring normal_cells may increase the correlation, tweaking dose and time does not. 
#ignoring neither blood cell or tnbc cells may increase the correlation
activity_RGES_summarized = aggregate(cbind(RGES, spearman, pearson, standard_value) ~ pert_iname, activity_RGES,  mean)
t.test(activity_RGES_summarized$RGES[activity_RGES_summarized$standard_value < 10000], activity_RGES_summarized$RGES[activity_RGES_summarized$standard_value > 10000])

spearman_test = cor.test(activity_RGES_summarized$spearman, log(activity_RGES_summarized$standard_value, 10), method="spearman")
pearson_test = cor.test(activity_RGES_summarized$pearson, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test = cor.test(activity_RGES_summarized$RGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
spearman_test
pearson_test
cor_test

lm_cmap_ic50 = lm(RGES ~ log(standard_value, 10), activity_RGES_summarized)
summary(lm_cmap_ic50)
summary(lm_cmap_ic50)$r.squared

pdf(paste( "fig/", cancer, "rges_ic50_min.pdf", sep=""))
ggplot(activity_RGES_summarized, aes(RGES, log(activity_RGES_summarized$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
        annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=2), sep=""), x = 0, y = 7.7, size = 6, colour = "black") +
        annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 7.3, size = 6, colour = "black") +
        scale_size(range = c(2, 5)) +
        xlab("RGES") + guides(shape=FALSE, size=FALSE) +
        ylab("log10(IC50) nm") + coord_cartesian(xlim = c(-1, 1), ylim=c(-1, 8)) 
dev.off()


###assess different RGES measures
activity_RGES_summarized = aggregate(cbind(RGES, spearman, pearson, standard_value) ~ pert_iname, activity_RGES,  min)
cor.test(activity_RGES_summarized$RGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor.test(activity_RGES_summarized$pearson, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor.test(activity_RGES_summarized$spearman, log(activity_RGES_summarized$standard_value, 10), method="spearman")

activity_RGES$activity = 1
activity_RGES$activity[activity_RGES$standard_value>10000] = 2

pred = prediction(activity_RGES$RGES, activity_RGES$activity)
performance(pred, "auc")
perf <- performance( pred, "tpr", "fpr" )
plot( perf )
abline(coef=c(0,1))

perf1 <- performance(pred, "prec", "rec")
plot(perf1)

rocobj <- roc(activity_RGES$activity,activity_RGES$RGES)
cutpoint <- coords(rocobj,x='local maximas',input='sensitivity', best.method = 'youden')
#cutpoint
#BRCA: 0.66 AUC
#threshold specificity sensitivity 
#-0.2729781   0.4461538   0.9318182 
#LIHC: 0.78
#-0.42 (local maximas)
#COAD 0.65
#-0.21


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

pdf(paste( "fig/", cancer, "rges_ic50_normalized.pdf", sep=""))
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

write.csv(activity_RGES_summarized, paste(cancer, "/rges_ic50_normalized.csv", sep=""))

pert_iname_resid = data.frame(resid = lm_cmap_ic50$residuals, activity_RGES_summarized)

pert_iname_resid = pert_iname_resid[order(abs(pert_iname_resid$resid)),]
tail(pert_iname_resid)

subset(lincs_drug_prediction, pert_iname == "vinblastine" & cell_id %in% cell_lines$LINCS)


########
#visualize drugs
tumor_drugs = read.csv("raw/tumor_drug_freq.csv")
tumor_drugs = tumor_drugs[tumor_drugs$cancer == cancer & tumor_drugs$Freq>5, ]

activity_RGES_summarized$clinic_used = match(activity_RGES_summarized$pert_iname, tumor_drugs$drugs)
activity_RGES_summarized$clinic_used = sapply(activity_RGES_summarized$clinic_used, function(x){ifelse(!is.na(x), "Therapy", "Others")})
activity_RGES_summarized$p_text = sapply(1:nrow(activity_RGES_summarized), function(id){
  if (activity_RGES_summarized$clinic_used[id] == "Therapy"){
    as.character(activity_RGES_summarized$pert_iname[id])
  }else{
    ""
  }
})
  
pdf(paste( "fig/", cancer, "rges_ic50_normalized_therapy.pdf", sep=""))
ggplot(activity_RGES_summarized, aes(sRGES, log(activity_RGES_summarized$standard_value, 10), label=p_text )) +  theme_bw()  +   geom_text(hjust=0, vjust=1, size=5) +
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18)) +                                                                                               
  stat_smooth(method="lm", se=F, color="black")  + geom_point(aes(color = clinic_used)) + 
  annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=2), sep=""), x = 0, y = 7.7, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 7.3, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("RGES") + guides(shape=FALSE, size=FALSE) +
  ylab("log10(IC50) nm") + coord_cartesian(xlim = c(-0.5, 0.5), ylim=c(-1, 8))
dev.off()

##check if correlation is retained without considering cell lines from the same lineage, 
#without same lineage
cell_lines = read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))
activity_RGES_subset = subset(activity_RGES, !(cell_id %in% cell_lines$LINCS)) 
results_subset = merge(activity_RGES_subset, cell_line_cancer, by.x="cell_id", by.y="lincs_cell_id")
results_subset$sRGES = sapply(1:nrow(results_subset), function(id){getsRGES3(results_subset$RGES[id], results_subset$cor[id], results_subset$pert_dose[id], results_subset$pert_time[id], diff, max(results_subset$cor))})

activity_RGES_summarized = aggregate(cbind(sRGES, standard_value) ~ pert_iname , results_subset, mean)
cor_test = cor.test(activity_RGES_summarized$sRGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test

#with the same lineage
cell_lines = read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))
activity_RGES_subset = subset(activity_RGES, (cell_id %in% cell_lines$LINCS)) 
results_subset = merge(activity_RGES_subset, cell_line_cancer, by.x="cell_id", by.y="lincs_cell_id")
results_subset$sRGES = sapply(1:nrow(results_subset), function(id){getsRGES3(results_subset$RGES[id], results_subset$cor[id], results_subset$pert_dose[id], results_subset$pert_time[id], diff, max(results_subset$cor))})

activity_RGES_summarized = aggregate(cbind(sRGES, standard_value) ~ pert_iname , results_subset, mean)
cor_test = cor.test(activity_RGES_summarized$sRGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test

