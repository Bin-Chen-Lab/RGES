# RGES and dose correlation
#some compounds may have been administrated longer or higher dose to exhibit effect

cancer = "BRCA"
lincs_drug_activity = read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)
lincs_drug_activity = aggregate(standard_value ~ pert_iname, lincs_drug_activity, median)
cell_lines = read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))

output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
lincs_drug_prediction = read.csv(output_path)

cdata <- ddply(lincs_drug_prediction, .(pert_iname, cell_id ), summarise, 
               time_count = length(unique(pert_time)),
               dose_count = length(unique(pert_dose))
)
cdata = cdata[order(cdata$dose_count), ]  

lincs_drug_prediction$RGES = lincs_drug_prediction$cmap_score
lincs_drug_prediction = subset(lincs_drug_prediction, cell_id == "MCF7" & pert_dose >0  & cell_id %in% cell_lines$LINCS)

#vorinostat geldanamycin trichostatin-a  gemcitabine  withaferin-a   wortmannin sirolimus wortmannin tretinoin tozasertib  
drugs = c("vorinostat", "geldanamycin", "trichostatin-a") #, "wortmannin", "tanespimycin", "withaferin-a") ,  "gemcitabine"

for (drug in drugs){
  for (time in c(6, 24)){
    lincs_drug_prediction_subset = subset(lincs_drug_prediction,  pert_iname == drug & pert_time == time)
    lincs_drug_prediction_subset = aggregate(RGES ~ pert_dose, lincs_drug_prediction_subset, mean)
    lincs_drug_prediction_subset = lincs_drug_prediction_subset[order(lincs_drug_prediction_subset$pert_dose),]
    cor_test = cor.test(lincs_drug_prediction_subset$RGES, log(lincs_drug_prediction_subset$pert_dose, 10), method="spearman")
    print(paste(drug, time, cor_test$p.value, cor_test$estimate))
    lm_dose_rges = lm(RGES ~ log(pert_dose, 10), lincs_drug_prediction_subset)
    
    pdf(paste( "fig/rges_dose_", cancer, drug, time, ".pdf", sep=""))
    print(ggplot(lincs_drug_prediction_subset, aes(y=RGES, x=log(pert_dose, 10)  )) +  theme_bw()  + 
            theme(legend.position ="bottom", axis.text=element_text(size=28), axis.title=element_text(size=28))  +                                                                                              
            stat_smooth(method="lm", se=F, color="black")  + geom_point() + 
            annotate("text", label = paste("r=", format(summary(lm_dose_rges)$r.squared ^ 0.5, digit=2), ", ",
                                           "P=", format(anova(lm_dose_rges)$`Pr(>F)`[1], digit=2), sep=""), 
                     x = 0, y = .20, size = 14, colour = "black") +
            annotate("text", label = paste(drug, ", ", time, "H", sep=""), x = 0, y = .3, size = 14, colour = "black") +
            scale_size(range = c(2, 5)) +
            xlab("log10(concentration) um") + guides(shape=FALSE, size=FALSE) +
            ylab("RGES") + coord_cartesian(xlim = c(-3, 2), ylim=c(-0.6, .4)) 
    )
    dev.off()
  }
}
