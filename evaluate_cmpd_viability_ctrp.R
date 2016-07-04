cancer = "BRCA"
site = "breast" #breast #large_intestine
drug_pred = read.csv(paste( cancer, "/lincs_all_score.csv", sep=""))
drug_pred$RGES = drug_pred$cmap_score
#drug_pred = subset(drug_pred, cell_id == "HEPG2")
load("~/Documents/stanford/sensitivity/data/ctrp/drug_cell_matrix.RData")

ctrp_cell_lines = read.csv("~/Documents/stanford/sensitivity/data/ctrp/cell_line.csv")
ctrp_cell_lines_site = ctrp_cell_lines$cell_line_name[ctrp_cell_lines$ccle_primary_site == site]

###
#evaluate shared cell lines individually
cell_lines = read.csv(paste("../table/", cancer, "_cell_lines.csv", sep=""))
cell_lines_subset = subset(cell_lines, LINCS != "" & CTRP !="")
cors = NULL
ps = NULL
for (i in 1:nrow(cell_lines_subset)){
  drug_pred_subset = subset(drug_pred, cell_id %in% cell_lines_subset$LINCS[i])
  drug_cell_matrix_subset =  drug_cell_matrix[,colnames(drug_cell_matrix) %in% cell_lines_subset$CTRP[i]]
  
  drug_auc = data.frame(drug = tolower(names(drug_cell_matrix_subset)), auc = drug_cell_matrix_subset)
  drug_auc = drug_auc[!is.na(drug_auc$auc),]
  
  drug_pred_subset = aggregate(RGES ~ pert_iname, drug_pred_subset, mean)
  drug_pred_subset$pert_iname = tolower(drug_pred_subset$pert_iname)
  activity_RGES_summarized = merge(drug_pred_subset, drug_auc, by.x="pert_iname", by.y="drug")
  test = cor.test(activity_RGES_summarized$RGES, activity_RGES_summarized$auc, method="spearman")
  plot(activity_RGES_summarized$RGES, activity_RGES_summarized$auc)
  cors = c(cors, test$estimate)
  ps = c(ps, test$p.value)
}

cell_cor = data.frame(cell_line = cell_lines_subset$LINCS, cor = cors, p = ps)

TNBC = c("MDAMB231", "HS578T", "SKBR3", "BT20")

drug_cell_matrix_subset =  drug_cell_matrix[,colnames(drug_cell_matrix) %in% cell_lines_subset$CTRP]
drug_auc = apply(drug_cell_matrix_subset, 1, function(values){median(values, na.rm=T)})
drug_auc = data.frame(drug = tolower(names(drug_auc)), auc = drug_auc)
drug_auc = drug_auc[!is.na(drug_auc$auc),]

drug_pred$pert_iname = tolower(drug_pred$pert_iname)

activity_RGES_summarized = merge(drug_pred, drug_auc, by.x="pert_iname", by.y="drug")
plot(activity_RGES_summarized$RGES, activity_RGES_summarized$auc)

cor.test(activity_RGES_summarized$RGES, activity_RGES_summarized$auc, method="spearman")

cor_test = cor.test(activity_RGES_summarized$RGES, activity_RGES_summarized$auc, method="spearman")
cor_test
lm_cmap_ic50 = lm(RGES ~ auc, activity_RGES_summarized)
summary(lm_cmap_ic50)
summary(lm_cmap_ic50)$r.squared

pdf(paste( "../fig/", cancer, "rges_auc_ctrp.pdf", sep=""))
p <- ggplot(activity_RGES_summarized, aes(RGES, activity_RGES_summarized$auc )) +  theme_bw()  + 
  theme(legend.position ="bottom", axis.text=element_text(size=18), axis.title=element_text(size=18))                                                                                                
print(p + stat_smooth(method="lm", se=F, color="black")  + geom_point() + 
        annotate("text", label = paste("r=", format(summary(lm_cmap_ic50)$r.squared ^ 0.5, digit=2), sep=""), x = 0, y = 18.5, size = 6, colour = "black") +
        annotate("text", label = paste("rho=", format(cor_test$estimate, digit=2), ", P=", format(cor_test$p.value, digit=3, scientific=T), sep=""), x = 0, y = 17, size = 6, colour = "black") +
        scale_size(range = c(2, 5)) +
        xlab("RGES") + guides(shape=FALSE, size=FALSE) +
        ylab("AUC") + coord_cartesian(xlim = c(-0.3, 0.3), ylim=c(0, 20))) 
dev.off()
