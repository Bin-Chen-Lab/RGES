########
##visualize IC50 for each cancer

cancer <- "LIHC"
lincs_drug_activity <- read.csv(paste(cancer, "/lincs_drug_activity_confirmed.csv", sep=""), stringsAsFactors=F)

drug_ic50_sd <- aggregate(standard_value ~ pref_name, lincs_drug_activity, sd)
drug_ic50_sd <- drug_ic50_sd[order(drug_ic50_sd$standard_value),]

pdf(paste(cancer,  "/lincs_drug_activity.pdf", sep=""))
par(mar = c(5, 4, 4, 2) + 0.5)

drug_order <- aggregate(standard_value ~ pref_name, lincs_drug_activity, median)
drug_ordered <- drug_order$pref_name[order(drug_order$standard_value, decreasing=T)]

lincs_drug_activity$pref_name <- factor(lincs_drug_activity$pref_name, levels = drug_ordered)

p <- ggplot(lincs_drug_activity, aes(((pref_name)), log(standard_value,10)))
print(p + geom_boxplot(notchwidth = 0.8) +  geom_jitter()  +  ylab("log(IC50), nm") +
        xlab("")  +
        geom_hline(yintercept = 4, colour="red", linetype = "longdash") +  theme_bw()  + 
        theme(axis.text.y = element_text(size=18),axis.title.y = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1, size=10))
)
dev.off()