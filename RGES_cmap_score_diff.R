###RGES vs conncectivity score
#
library(ggplot2)
#RGES score distribution
output_path <- paste("BRCA", "/all_lincs_score.csv", sep="")
lincs_drug_prediction_new <- read.csv(output_path)
lincs_drug_prediction_new <- subset(lincs_drug_prediction_new, pert_dose >0)
lincs_drug_prediction_new$RGES <- lincs_drug_prediction_new$cmap_score
lincs_drug_prediction_new$type <- "RGES"

output_path <- paste("BRCA", "/all_lincs_score_cmap_old.csv", sep="")
lincs_drug_prediction_old <- read.csv(output_path)
lincs_drug_prediction_old <- subset(lincs_drug_prediction_old, pert_dose >0)
lincs_drug_prediction_old$RGES <- lincs_drug_prediction_old$cmap_score
lincs_drug_prediction_old$type <- "Connectivity Score"

lincs_drug_predictions <- rbind(lincs_drug_prediction_new, lincs_drug_prediction_old)
pdf("fig/rges_cmap_score_diff.pdf")
ggplot(lincs_drug_predictions, aes(x=RGES)) + geom_density(aes(group=type, fill=type, colour=type),alpha = 0.5) + 
  theme_bw() +  
  theme(text = element_text(size=24), 
        panel.background = element_rect(color = 'white'), axis.text = element_text( size=24), 
        axis.title.x = element_text( size=24), axis.title.y = element_text( size=24),
        legend.position=c(0.2, 0.9), legend.title =  element_blank()
  ) 
dev.off()


#correlation with IC50
lincs_drug_activity <- read.csv("BRCA/rges_ic50_normalized.csv")
lincs_drug_prediction_new_ic50 <- merge(lincs_drug_prediction_new, lincs_drug_activity, by="pert_iname")
lincs_drug_prediction_old_ic50 <- merge(lincs_drug_prediction_old, lincs_drug_activity, by="pert_iname")


lincs_drug_prediction_new_ic50_merged <- aggregate(cbind(cmap_score, standard_value) ~ pert_iname,lincs_drug_prediction_new_ic50, median )
pdf("fig/rges_cmap_score_diff_cor_ic50_cmap.pdf")
ggplot(lincs_drug_prediction_new_ic50_merged, aes(cmap_score, log(lincs_drug_prediction_new_ic50_merged$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18)) +      
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  xlab("RGES") +
  ylab("log10(IC50) nm") 
dev.off()

lincs_drug_prediction_old_ic50_merged <- aggregate(cbind(cmap_score, standard_value) ~ pert_iname,lincs_drug_prediction_old_ic50, median )
pdf("fig/rges_cmap_score_diff_cor_ic50_rges.pdf")
ggplot(lincs_drug_prediction_old_ic50_merged, aes(cmap_score, log(lincs_drug_prediction_old_ic50_merged$standard_value, 10)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18)) +      
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  xlab("RGES") +
  ylab("log10(IC50) nm") 
dev.off()

