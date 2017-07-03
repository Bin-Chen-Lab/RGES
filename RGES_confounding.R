#RGES depends on cell line, dose, and time

library(plyr)
library(ggplot2)

cancers <- c("BRCA", "LIHC", "COAD") #, "random"
landmark <- 1

#RGES distribution in top commonly profiled compounds
for (cancer in cancers[1:3]){
  output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
  lincs_drug_prediction <- read.csv(output_path)
  lincs_drug_prediction <- subset(lincs_drug_prediction, pert_dose >0)
  
  cell_lines <- read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))
  
  lincs_drug_prediction$same_lineage <- F
  lincs_drug_prediction$same_lineage[lincs_drug_prediction$cell_id  %in% cell_lines$LINCS] <- T
  
  lincs_drug_prediction$RGES <- lincs_drug_prediction$cmap_score
  cdata_all <- ddply(lincs_drug_prediction, .(pert_iname ), summarise, 
                 n = length(RGES),
                 min = min(RGES), 
                 max = max(RGES)
  )
  mean(cdata_all$n)
  median(cdata_all$n)
  sd(cdata_all$n)
  
  cdata <- cdata_all[cdata_all$n>1, ]
  cdata <- cdata[order(cdata$n, decreasing=T),]
  cdata <- cdata[1:4, ]
  cdata <- cdata[order(cdata$n, decreasing=F),]
  
  lincs_drug_prediction_subset <- merge(lincs_drug_prediction, cdata, by="pert_iname")
  lincs_drug_prediction_subset$pert_iname <- factor(lincs_drug_prediction_subset$pert_iname, levels = cdata$pert_iname)
 
  pdf(paste("fig/", cancer, "rges_top5.pdf", sep=""))
  
  print(
   ggplot(lincs_drug_prediction_subset, aes(x = pert_iname, 
    y = RGES, colour = same_lineage)) + geom_point(size = 0.5, position = position_jitter(width = 0.5, 
        height = 0)) + 
    theme_bw() + coord_flip() + ylim(c(-0.8,0.8)) + 
    theme(legend.position = "bottom", text = element_text(size=18), axis.text = element_text( size=18))  +
     scale_color_manual(name = "", values = c("TRUE" = "red", "FALSE" = "gray"), labels = c( "others", paste(cancer, "cell lines"))) + 
    theme(legend.title = element_text(size = 18), plot.margin=unit(c(9,1,1,1), "cm")) +  xlab("") + ylab("") + 
     annotate("text", label = paste(cancer), 
              x = 1, y = 0.7, size = 6, colour = "black") 
  )
  dev.off()

} 


#####
###
lincs_drug_predictions <- data.frame()
for (cancer in cancers){
  output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
  lincs_drug_prediction <- read.csv(output_path)
  lincs_drug_prediction <- subset(lincs_drug_prediction, pert_dose >0)
  
  lincs_drug_prediction$RGES <- lincs_drug_prediction$cmap_score
  lincs_drug_prediction$cancer <- cancer
  lincs_drug_predictions <- rbind(lincs_drug_predictions, lincs_drug_prediction)
}
lincs_drug_predictions$cell_type <- paste(lincs_drug_predictions$cancer, "all cell lines", sep=", ")

#from same lineage (random gene set has no cell line)
lincs_drug_predictions_same_lineage <- data.frame()
for (cancer in cancers[1:3]){
  output_path <- paste(cancer, "/all_lincs_score.csv", sep="")
  lincs_drug_prediction <- read.csv(output_path)
  lincs_drug_prediction <- subset(lincs_drug_prediction, pert_dose >0)
  
  cell_lines <- read.csv(paste("raw/cell_lines/", cancer, "_cell_lines.csv", sep=""))
  lincs_drug_prediction <- subset(lincs_drug_prediction, cell_id %in% cell_lines$LINCS)
  lincs_drug_prediction$RGES <- lincs_drug_prediction$cmap_score
  lincs_drug_prediction$cancer <- cancer
  lincs_drug_predictions_same_lineage <- rbind(lincs_drug_predictions_same_lineage, lincs_drug_prediction)
}
lincs_drug_predictions_same_lineage$cell_type <- paste(lincs_drug_predictions_same_lineage$cancer, "same lineage", sep=", ")

pdf("fig/rges_dist.pdf")
ggplot(lincs_drug_predictions, aes(x=RGES)) + geom_density(aes(group=cancer, fill=cancer, colour=cancer),alpha = 0.5) + 
  theme_bw() +  
  theme(text = element_text(size=18), 
        panel.background = element_rect(color = 'white'), axis.text = element_text( size=18), 
        axis.title.x = element_text( size=20), axis.title.y = element_text( size=20),
        legend.position=c(0.8, 0.8), legend.title =  element_blank()
  ) 
dev.off()

pdf("fig/rges_dist_same_lineage.pdf")
ggplot(rbind(lincs_drug_predictions_same_lineage), aes(x=RGES)) + geom_density(aes(group=cell_type, fill=cancer, colour=cancer),alpha = 0.5) + 
  theme_bw() +  
  theme(text = element_text(size=18), 
        panel.background = element_rect(color = 'white'), axis.text = element_text( size=18), 
        axis.title.x = element_text( size=20), axis.title.y = element_text( size=20),
        legend.position=c(0.8, 0.8), legend.title =  element_blank()
  ) 
dev.off()

#variation of score across cell lines
cdata <- ddply(lincs_drug_predictions_same_lineage, .(pert_iname, pert_dose, pert_time, cancer ), summarise, 
               n = length(RGES),
               mean = mean(RGES),
               sd = sd(RGES)
)
cdata <- subset(cdata, !is.na(sd))        

tapply(cdata$sd, cdata$cancer, summary)


#sd accross multiple cell lines from the same lineage
pdf("fig/rges_sd_cells_same_lineage.pdf")
ggplot(cdata, aes(x=sd)) + geom_density(aes(group=cancer, fill=cancer, colour=cancer), alpha=0.5) + theme_bw() + xlab("s.d. of RGES of individual compounds across cell lines") +
  theme(text = element_text(size=18), 
        panel.background = element_rect(color = 'white'), axis.text = element_text( size=18), 
        axis.title.x = element_text( size=20), axis.title.y = element_text( size=20),
        legend.position=c(0.8, 0.8), legend.title =  element_blank()
  ) 
dev.off()

hist(cdata$sd)             
summary(cdata$sd)

mean(lincs_drug_prediction$RGES)

#variation of score across cell lines
cdata <- ddply(lincs_drug_predictions, .(pert_iname, pert_dose, pert_time,  cancer ), summarise, 
               n = length(RGES),
               mean = mean(RGES),
               sd = sd(RGES)
)
cdata = subset(cdata, !is.na(sd))        

tapply(cdata$sd, cdata$cancer, summary)

#sd accross multiple cells 
pdf("fig/rges_sd_cells.pdf")
ggplot(cdata, aes(x=sd)) + geom_density(aes(group=cancer, fill=cancer, colour=cancer), alpha=0.5) + theme_bw() + xlab("s.d. of RGES of individual compounds across cell lines") +
  theme(text = element_text(size=18), 
        panel.background = element_rect(color = 'white'), axis.text = element_text( size=18), 
        axis.title.x = element_text( size=20), axis.title.y = element_text( size=20),
        legend.position=c(0.8, 0.8), legend.title =  element_blank()
  ) 
dev.off()

#cross plate or replicates
cdata_plate <- ddply(lincs_drug_predictions, .(pert_iname, pert_dose, pert_time, cell_id, cancer ), summarise, 
               n = length(RGES),
               mean = mean(RGES),
               sd = sd(RGES)
)
cdata_plate <- subset(cdata_plate, !is.na(sd))        

cdata$cell_id <- NA
cdata_plate$type <- "same"
cdata$type <- "different"
cdata_plate_cell <- rbind(cdata, cdata_plate)

cdata_plate_cell$cancer <- factor(cdata_plate_cell$cancer, levels = c("BRCA", "LIHC", "COAD"))

pdf("fig/rges_cell.pdf")
ggplot(cdata_plate_cell, aes(y=sd, x=factor(cancer), fill=factor(type))) + geom_boxplot(alpha=0.5)+ ylab("s.d. of RGES within same compounds") + xlab("") + theme_bw()  +
  theme(text = element_text(size=18), 
        panel.background = element_rect(color = 'white'), axis.text = element_text( size=18), 
        axis.title.x = element_text( size=20), axis.title.y = element_text( size=20),
        legend.position="top") +   
  scale_fill_manual(name="cell line", values = c("red", "black")) 
dev.off()


#same drug same cell line and dose, different treatment duration
drug_6_24s <- data.frame()
comparison_time_results <- NULL
for (cancer_type in cancers){
  lincs_drug_prediction_6h <- subset(lincs_drug_predictions, pert_time < 24 & cancer == cancer_type)
  lincs_drug_prediction_24h <- subset(lincs_drug_predictions, pert_time >= 24 & cancer == cancer_type)

  drug_6_24 <- merge(lincs_drug_prediction_6h, lincs_drug_prediction_24h, by=c("pert_iname", "cell_id", "pert_dose"))
  
  lm_model <- lm(RGES.y ~ RGES.x + pert_dose, drug_6_24)
  
  drug_6_24s <- rbind(drug_6_24s, data.frame(RGES = drug_6_24$RGES.x, seq=seq(1, nrow(drug_6_24)), cancer = cancer_type, time = "<24h"))
  drug_6_24s <- rbind(drug_6_24s, data.frame(RGES = drug_6_24$RGES.y, seq=seq(1, nrow(drug_6_24)), cancer = cancer_type, time = ">=24h"))
    
  test <- t.test(drug_6_24$RGES.x, drug_6_24$RGES.y, paired =T)$p.value
  ratio <- mean(drug_6_24$RGES.y - drug_6_24$RGES.x)
  longer <- sum((drug_6_24$RGES.y - drug_6_24$RGES.x)<0)
  shorter <- sum((drug_6_24$RGES.y - drug_6_24$RGES.x)>0)
  comparison_time_results <- c(comparison_time_results, list(test, ratio, longer, shorter))
}

pdf("fig/rges_time.pdf")
ggplot(drug_6_24s, aes(y=RGES, x=factor(cancer), fill=factor(time))) + geom_boxplot(alpha=0.5)+ xlab("") + theme_bw()  +
  theme(text = element_text(size=18), 
        panel.background = element_rect(color = 'white'), axis.text = element_text( size=18), 
        axis.title.x = element_text( size=20), axis.title.y = element_text( size=20),
        legend.position="top") +   
        scale_fill_manual(name="Treatment duration", values = c("green", "pink")) 
dev.off()

drug_6_24s_subset <- subset(drug_6_24s, cancer == "LIHC")

sum(drug_6_24s_subset$RGES.x > drug_6_24s_subset$RGES.y)
sum(drug_6_24s_subset$RGES.x < drug_6_24s_subset$RGES.y)

summary(drug_6_24$RGES.x- drug_6_24$RGES.y)

#drug_6_24[(drug_6_24$RGES.x - drug_6_24$RGES.y) < -0.4, ]

#dose 
drug_dose_1_2s <- data.frame()
comparison_dose_results <- NULL
for (cancer_type in cancers){
  lincs_drug_prediction_dose1 <- subset(lincs_drug_predictions, pert_dose < 10   & cancer == cancer_type)
  lincs_drug_prediction_dose2 <- subset(lincs_drug_predictions, pert_dose >= 10  & cancer == cancer_type)
  
  drug_dose_1_2 <- merge(lincs_drug_prediction_dose1, lincs_drug_prediction_dose2, by=c("pert_iname", "cell_id", "pert_time"))
  
  #lm_model <- glm(RGES.y ~ RGES.x + pert_time, family = "gaussian", drug_dose_1_2)
  
  drug_dose_1_2s <- rbind(drug_dose_1_2s, data.frame(RGES = drug_dose_1_2$RGES.x, cancer = cancer_type, dose = "<10um"))
  drug_dose_1_2s <- rbind(drug_dose_1_2s, data.frame(RGES = drug_dose_1_2$RGES.y, cancer = cancer_type, dose = ">=10um"))
  
  test <- t.test(drug_dose_1_2$RGES.x, drug_dose_1_2$RGES.y, paired =T)$p.value
  ratio <- mean(drug_dose_1_2$RGES.y- drug_dose_1_2$RGES.x)
  comparison_dose_results <- c(comparison_dose_results, list(test, ratio))
}

pdf("fig/rges_dose.pdf")
ggplot(drug_dose_1_2s, aes(y=RGES, x=factor(cancer), fill=factor(dose))) + geom_boxplot(alpha=0.5) + xlab("") + theme_bw()  +
  theme(text = element_text(size=18), 
        panel.background = element_rect(color = 'white'), axis.text = element_text( size=18), 
        axis.title.x = element_text( size=20), axis.title.y = element_text( size=20),
        legend.position="top") +  
        scale_fill_manual(name="Treatment concentration", values = c("blue", "yellow")) 
dev.off()

#
drug_class <- read.csv("~/Documents/stanford/sensitivity/data/GDSC1000/compounds.csv", stringsAsFactors = F)
cytotoxic <- drug_class$Name[drug_class$Action == "cytotoxic"]
targeted <- drug_class$Name[drug_class$Action == "targeted"]

#identify drugs whose RGES were not affected by the confounding factors?
invariant_cmpds <- NULL
variant_cmpds <- NULL
for (cancer_type in cancers){
  #cancer_type <- "BRCA"
  lincs_drug_prediction_subset <- subset(lincs_drug_predictions,  cancer == cancer_type)
  cmpd_rges_variation <- by(lincs_drug_prediction_subset$RGES, lincs_drug_prediction_subset$pert_iname, sd)
  cmpd_rges_dose <- by(lincs_drug_prediction_subset$pert_dose, lincs_drug_prediction_subset$pert_iname, function(x) length(unique(x)))
  cmpd_rges_time <- by(lincs_drug_prediction_subset$pert_time, lincs_drug_prediction_subset$pert_iname, function(x) length(unique(x)))
  cmpd_rges_cell_id <- by(lincs_drug_prediction_subset$cell_id, lincs_drug_prediction_subset$pert_iname, function(x) length(unique(x)))
  cmpd_rges_count <- by(lincs_drug_prediction_subset$RGES, lincs_drug_prediction_subset$pert_iname, length)
  
  cmpd_rges_variation_subset <- cmpd_rges_variation[cmpd_rges_dose>1 & cmpd_rges_time > 1 & cmpd_rges_cell_id>1]
  cmpd_rges_variation_subset <- cmpd_rges_variation_subset[!is.na(cmpd_rges_variation_subset)]
  
  head(sort(cmpd_rges_variation_subset), 30)
  tail(sort(cmpd_rges_variation_subset), 30)
  cutoff1 <- qnorm(0.01, mean(cmpd_rges_variation_subset), sd(cmpd_rges_variation_subset))
  cutoff2 <- qnorm(0.01, mean(cmpd_rges_variation_subset), sd(cmpd_rges_variation_subset), lower.tail = F)
  
  invariant_cmpds <- c(invariant_cmpds, names(cmpd_rges_variation_subset[cmpd_rges_variation_subset < cutoff1]))
  variant_cmpds <- c(variant_cmpds, names(cmpd_rges_variation_subset[cmpd_rges_variation_subset > cutoff2]))
  
  cmpd_rges_variation_cytotoxic <- cmpd_rges_variation_subset[tolower(names(cmpd_rges_variation_subset)) %in% tolower(cytotoxic)]
  cmpd_rges_variation_targeted <- cmpd_rges_variation_subset[tolower(names(cmpd_rges_variation_subset)) %in% tolower(targeted)]
  print(t.test(cmpd_rges_variation_targeted, cmpd_rges_variation_cytotoxic, alternative = "greater"))
}

sort(table(invariant_cmpds))
sort(table(variant_cmpds))

qnorm(0.01, mean(cmpd_rges_variation_subset), sd(cmpd_rges_variation_subset), lower.tail = F)

############
#explore machine learning methods to adjust RGES
###########
#find drugs tested at least in two distinct dose, two distinct times, and two distinct cell lines
cmpds <- names(cmpd_rges_variation_subset)
p <- sapply(cmpds, function(cmpd){
  a <- lincs_drug_prediction_subset[lincs_drug_prediction_subset$pert_iname == cmpd, ]
  rges_lm <- lm(RGES ~ pert_dose + pert_time + cell_id , a)
  #rges_lm_summary <- summary(rges_lm)
  #anova(rges_lm)$`Pr(>F)`[1]
  f <- summary(rges_lm)$fstatistic
  pf(f[1],f[2],f[3],lower.tail=F)
})

cmpd_p <- data.frame(cmpd = cmpds, p)
cmpd_p <- cmpd_p[order(cmpd_p$p), ]
sum(cmpd_p$p < 0.05, na.rm = T)
sum(cmpd_p$p > 0.05, na.rm = T)

a <- lincs_drug_prediction_subset[lincs_drug_prediction_subset$pert_iname == "oxetane", ]
rges_lm <- lm(RGES ~  pert_dose +  pert_time + cell_id , a)
summary(rges_lm)
#anova(rges_lm)$`Pr(>F)`
#names((rges_lm)) 

rges_glm <- glm(RGES ~ pert_dose + pert_time + cell_id , family = "gaussian", a)
summary(rges_glm)

library("Rcmdr")
leveneTest(a$RGES, lincs_drug_prediction_subset$RGES)

