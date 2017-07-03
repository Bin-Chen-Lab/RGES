#compute correlation between RGES and drug IC50s in single cell lines
#using drug IC50 from external databases (CCLE, GDSC)
library(cowplot)

######################
#MAIN
#####################
###
cancer <- "BRCA"
cell_line_selected <- "MCF7" #HT29, HEPG2, MCF7
landmark <- 1

#build a inference model according to dose and time
output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
lincs_drug_prediction <- read.csv(output_path)

lincs_drug_prediction <- subset(lincs_drug_prediction,  pert_dose > 0 & pert_time %in% c(6, 24))
lincs_drug_prediction$RGES <- lincs_drug_prediction$cmap_score

lincs_drug_prediction_subset <- subset(lincs_drug_prediction, cell_id %in% c(cell_line_selected)) #HT29 MCF7
#lincs_drug_prediction_subset <- aggregate(RGES ~ pert_iname, lincs_drug_prediction_subset, median)
#by default, use this one
lincs_drug_prediction_subset$sRGES <- lincs_drug_prediction_subset$RGES

#################
###sensitivity data from CTRP
load("~/Documents/stanford/sensitivity/data/ctrp/drug_cell_matrix.RData")

ctrp_drug_activity_subset <- drug_cell_matrix[, cell_line_selected ]
ctrp_drug_activity_subset <- data.frame(pert_iname_lower = tolower(names(ctrp_drug_activity_subset)), standard_value = as.numeric(ctrp_drug_activity_subset))
ctrp_drug_activity_subset <- subset(ctrp_drug_activity_subset, !is.na(standard_value))
lincs_drug_prediction_subset$pert_iname_lower <- tolower(lincs_drug_prediction_subset$pert_iname)

activity_RGES <- merge(lincs_drug_prediction_subset, ctrp_drug_activity_subset, by="pert_iname_lower")

activity_RGES_summarized <- aggregate(cbind(sRGES, standard_value) ~ pert_iname, activity_RGES,  median)

activity_RGES_summarized$activity <- "effective"
activity_RGES_summarized$activity[activity_RGES_summarized$standard_value>10] <- "ineffective"
efficacy_test <- t.test(activity_RGES_summarized$sRGES[activity_RGES_summarized$activity == "effective"], activity_RGES_summarized$sRGES[activity_RGES_summarized$activity == "ineffective"])

cor_test <- cor.test(activity_RGES_summarized$sRGES, (activity_RGES_summarized$standard_value), method="spearman")
cor_test
lm_cmap_ic50 <- lm( standard_value ~ sRGES, activity_RGES_summarized)
summary(lm_cmap_ic50)
summary(lm_cmap_ic50)$r.squared^0.5

pdf(paste( "fig/", cancer, "rges_ic50_", cell_line_selected, "_normalized_CTRP.pdf", sep=""))
lm_plot <- ggplot(activity_RGES_summarized, aes(sRGES, (activity_RGES_summarized$standard_value)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label = paste(cancer, ",", cell_line_selected, sep=""), 
           x = 0, y = 18.1 + 2, size = 6, colour = "black") +
  annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=2), sep=""), 
           x = 0, y = 17.4 + 2, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 16.5 + 2, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("RGES") + guides(shape=FALSE, size=FALSE) +
  ylab("AUC") + coord_cartesian(xlim = c(-1, 1), ylim=c(5, 20)) 

bar_plot <-  ggplot(activity_RGES_summarized, aes(activity, sRGES  )) + geom_boxplot() +  coord_flip()

ggdraw() +
  draw_plot(lm_plot) 

dev.off()

########
#GDSC database
gdsc_drug_activity <- read.csv("raw/external_sensitivity/cmpd_ic50.csv", stringsAsFactors=F, row.names = 1)
gdsc_drug_activity_subset <- lincs_drug_activity[lincs_drug_activity$cell_line == cell_line_selected, ]
gdsc_drug_activity_subset <- data.frame(pert_iname_lower = tolower(colnames(gdsc_drug_activity_subset)), standard_value = as.numeric(gdsc_drug_activity_subset))
gdsc_drug_activity_subset <- subset(gdsc_drug_activity_subset, !is.na(standard_value))
gdsc_drug_activity_subset$pert_iname_lower <- sapply(gdsc_drug_activity_subset$pert_iname_lower, function(x) paste(unlist(strsplit(x, "\\.")), collapse = " "))
lincs_drug_prediction_subset$pert_iname_lower <- tolower(lincs_drug_prediction_subset$pert_iname)

activity_RGES <- merge(lincs_drug_prediction_subset, gdsc_drug_activity_subset, by="pert_iname_lower")

activity_RGES_summarized <- aggregate(cbind(sRGES, standard_value) ~ pert_iname, activity_RGES,  mean)

activity_RGES_summarized$activity <- "effective"
activity_RGES_summarized$activity[activity_RGES_summarized$standard_value>10] <- "ineffective"
efficacy_test <- t.test(activity_RGES_summarized$sRGES[activity_RGES_summarized$activity == "effective"], activity_RGES_summarized$sRGES[activity_RGES_summarized$activity == "ineffective"])

cor_test <- cor.test(activity_RGES_summarized$sRGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
cor_test
lm_cmap_ic50 <- lm( standard_value ~ sRGES, activity_RGES_summarized)
summary(lm_cmap_ic50)
summary(lm_cmap_ic50)$r.squared^0.5

pdf(paste( "fig/", cancer, "rges_ic50_", cell_line_selected, "_normalized_GDSC.pdf", sep=""))
lm_plot <- ggplot(activity_RGES_summarized, aes(sRGES, (activity_RGES_summarized$standard_value)  )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))  +                                                                                              
  stat_smooth(method="lm", se=F, color="black")  + geom_point(size=3) + 
  annotate("text", label = paste(cancer, ",", cell_line_selected, sep=""), 
           x = 0, y = 8.1, size = 6, colour = "black") +
  annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), ", ",  "P=", format(anova(lm_cmap_ic50)$`Pr(>F)`[1], digit=2), sep=""), 
           x = 0, y = 7.7, size = 6, colour = "black") +
  annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 7.3, size = 6, colour = "black") +
  scale_size(range = c(2, 5)) +
  xlab("RGES") + guides(shape=FALSE, size=FALSE) +
  ylab("IC50 (um)") + coord_cartesian(xlim = c(-1, 1), ylim=c(-1, 8)) 

bar_plot <-  ggplot(activity_RGES_summarized, aes(activity, sRGES  )) + geom_boxplot() +  coord_flip()

ggdraw() +
  draw_plot(lm_plot)

dev.off()

ctrp_gdsc <- merge(gdsc_drug_activity_subset, ctrp_drug_activity_subset, by = "pert_iname_lower")
cor.test(ctrp_gdsc$standard_value.x, ctrp_gdsc$standard_value.y)
plot(ctrp_gdsc$standard_value.x, ctrp_gdsc$standard_value.y)
