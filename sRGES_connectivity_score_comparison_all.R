brca = read.csv( paste("BRCA", "/compare_sRGES_lincs_cloud_", "reduced", ".csv", sep=""))
lihc = read.csv( paste("LIHC", "/compare_sRGES_lincs_cloud_", "reduced", ".csv", sep=""))
coad = read.csv( paste("COAD", "/compare_sRGES_lincs_cloud_", "reduced", ".csv", sep=""))

all <- cbind(brca[,2], lihc[,2], coad[,2])

colnames(all) <- c("BRCA", "LIHC", "COAD")
rownames(all) <- brca$X

library(pheatmap)

my_palette <- colorRampPalette(c("green", "white", "red"))(n = 50)
pheatmap(all, color = my_palette, cluster_rows=F, cluster_cols=F, 
         cellwidth=25,  cellheight=25, display_numbers=T, file=paste("fig/", "srges_lincs_cloud_performance.pdf", sep=""))
