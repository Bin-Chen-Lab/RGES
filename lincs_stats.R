#lincs cmpd/

library ("ggplot2")

cancer <- "BRCA"

#cmpd
output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
lincs_drug_prediction <- read.csv(output_path)

#~70 signatures have no dose info. Ingore them in the following analysis
lincs_drug_prediction <- subset(lincs_drug_prediction, pert_dose >0)

nrow(lincs_drug_prediction)
length(unique(lincs_drug_prediction$pert_iname))

table(lincs_drug_prediction$pert_time)

sum(lincs_drug_prediction$pert_dose == 10)
sum(lincs_drug_prediction$pert_dose > 1 & lincs_drug_prediction$pert_dose<10)
sum(lincs_drug_prediction$pert_dose < 1 & lincs_drug_prediction$pert_dose<10)

sum(lincs_drug_prediction$padj<0.05)/nrow(lincs_drug_prediction)
length(unique(lincs_drug_prediction$pert_iname[lincs_drug_prediction$padj<0.05]))/length(unique(lincs_drug_prediction$pert_iname))

lincs_drug_prediction$RGES <- lincs_drug_prediction$cmap_score
lincs_drug_prediction$dose[lincs_drug_prediction$pert_dose >=10] <- ">=10"
lincs_drug_prediction$dose[lincs_drug_prediction$pert_dose <10] <- "<10"

#visualize the top 10 cells
cells <- tail(sort(table(lincs_drug_prediction$cell_id)), 15)
lincs_drug_prediction$cell_line <- as.character(lincs_drug_prediction$cell_id)
lincs_drug_prediction$cell_line[!(lincs_drug_prediction$cell_line %in% names(cells))] <- "other 56 cell lines"
lincs_drug_prediction$cell_line <- factor(lincs_drug_prediction$cell_line, levels = c(rev(names(cells)), "other 56 cell lines"))

#
drug_freq <- data.frame((table(lincs_drug_prediction$pert_iname)))
#pert_iname distribution
pdf(paste( "fig/",  "pert_iname.pdf", sep=""))
ggplot(drug_freq, aes( Freq, fill="red", color="red")) + geom_histogram(bins=20) + theme_bw() +  scale_x_log10() + 
    theme(text = element_text(size = 18), 
          panel.background = element_rect(color = "white"), 
          axis.text = element_text(size = 18), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20), 
          legend.position = "none", legend.title = element_blank()) + 
  theme(panel.background = element_rect(fill = NA), 
        plot.background = element_rect(linetype = "blank")) + 
  labs(x = "# profiles per compound", y= "count") + ylim(0, 5500)
dev.off()

#dose distribution
pdf(paste( "fig/",  "dose.pdf", sep=""))
ggplot(lincs_drug_prediction, aes(x = dose, 
                                  fill = cell_line)) + geom_bar() + theme_bw() + 
  xlab("treatment concentration") + theme(text = element_text(size = 18), 
                                         panel.background = element_rect(color = "white"), 
                                         axis.text = element_text(size = 18), 
                                         axis.title.x = element_text(size = 20), 
                                         axis.title.y = element_text(size = 20), 
                                         legend.position = "none", legend.title = element_blank()) + 
  theme(panel.background = element_rect(fill = NA), 
        plot.background = element_rect(linetype = "blank")) + 
  labs(x = "treatment concentration (um)", y= "# drug expression profiles") + ylim(0, 50000)
dev.off()

pdf(paste( "fig/",  "pert_dose.pdf", sep=""))
ggplot(lincs_drug_prediction, aes(x = pert_dose
                                 )) + geom_histogram(binwidth=5) + theme_bw() + 
  xlab("treatment concentration") + theme(text = element_text(size = 18), 
                                          panel.background = element_rect(color = "white"), 
                                          axis.text = element_text(size = 18), 
                                          axis.title.x = element_text(size = 20), 
                                          axis.title.y = element_text(size = 20), 
                                          legend.position = "none", legend.title = element_blank()) + 
  theme(panel.background = element_rect(fill = NA), 
        plot.background = element_rect(linetype = "blank")) + 
  labs(x = "treatment concentration (um)", y= "# drug expression profiles") + ylim(0, 50000)
dev.off()

#time distribution
pdf(paste( "fig/",  "time.pdf", sep=""))
ggplot(lincs_drug_prediction, aes(x = as.factor(pert_time), 
    fill = cell_line)) + geom_bar() + theme_bw() + 
    xlab("treatment duration (H)") + theme(text = element_text(size = 18), 
    panel.background = element_rect(color = "white"), 
    axis.text = element_text(size = 18), 
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20), 
    legend.position = c(0.8, 0.6), legend.title = element_blank()) + 
    theme(panel.background = element_rect(fill = NA), 
        plot.background = element_rect(linetype = "blank")) + 
    labs(x = "treatment duration (h)", y= "# drug expression profiles") + ylim(0, 50000)

dev.off()

#cell line
pdf(paste( "fig/",  "cell_line.pdf", sep=""))
ggplot(lincs_drug_prediction, aes(x=as.factor(cell_line))) + geom_bar(fill="red", color="red") + theme_light() + xlab("") + 
  theme(text = element_text(size=18), 
        panel.background = element_rect(color = 'white'), axis.text = element_text( size=18), 
        axis.title.x = element_text( size=20), axis.title.y = element_text( size=20),
        legend.position=c(0.2, 0.8), legend.title =  element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
  ) +    labs( y= "# drug expression profiles") + ylim(0, 50000)

dev.off()


