#need to download data from GSE69845

# deal negative values
library(preprocessCore)

expr = read.delim("GEO/GSE69845_series_matrix.txt", skip = 68, header=T, sep="\t") #MCF7: GSE69845 #HEPG: GSE69850
label = read.csv("GEO/GSE69845_outcome.csv")

#only keep dose 10um and 
#hepg2_label = subset(hepg2_label, dose %in% c("10uM", "DMSO"))


  #library("GEOquery")
  gpl <- read.delim("GEO/GPL13667-15572.txt", sep="\t", skip = 43, header=T)
  geneMappings <- gpl[,c("ID","Entrez.Gene")]
  names(geneMappings) <- c("probe","GeneID")
  geneMappings <- subset(geneMappings, !is.na(GeneID))
  dupProbes <- unique(geneMappings$probe[duplicated(geneMappings$probe)])
  geneMappings <- geneMappings[! geneMappings$probe %in% dupProbes,] #remove duplicated probes
  rownames(geneMappings) <- geneMappings$probe

  #remove probes if half of values are missed
  expr_gene = merge(geneMappings, expr,  by.x="probe", by.y="ID_REF")
  probes = expr_gene$GeneID
  comparison_frame = as.matrix(expr[, -c(1,2)])
  complete_probes <- apply(comparison_frame,1,function(row) {
    ifelse (sum(is.na(row))>10, FALSE, TRUE)
  })
  comparison_frame <- comparison_frame[complete_probes,]
  probes <- probes[complete_probes ]
  

  comparison_frame <- log2(comparison_frame+0.001) # add 1 to avoid log2(0)

  #normalization
  #boxplot(comparison_frame)
  comparison_frame1 <- normalize.quantiles(as.matrix(comparison_frame ))
  row.names(comparison_frame1) <- row.names(comparison_frame)
  colnames(comparison_frame1) <- colnames(comparison_frame)
  comparison_frame <- comparison_frame1
  
  comparison_frame_merged <- aggregate(. ~ probes, data.frame(probes,comparison_frame ), mean)

  drugs = as.character(unique(label$drug))
  controls = as.character(label$GSM[label$drug == "DMSO"])
  
  #cor(comparison_frame_merged[, colnames(comparison_frame_merged) %in% controls])
  control_expr = apply(comparison_frame_merged[, colnames(comparison_frame_merged) %in% controls], 1, mean)
  drug_sig = data.frame(GeneID=comparison_frame_merged$probes)
  for (drug in drugs){
    treatment_expr = apply(comparison_frame_merged[, colnames(comparison_frame_merged) %in% as.character(label$GSM[label$drug == drug])], 1, mean)
    drug_sig = cbind(drug_sig, treatment_expr - control_expr)
  }
  colnames(drug_sig) = c("GeneID", drugs)
  drug_sig = drug_sig[drug_sig$GeneID != "---", ]
  drug_sig$GeneID = sapply(as.character(drug_sig$GeneID), function(id){
    unlist(strsplit(id, " /// "))[1]
  })
  
  drug_sig$GeneID = as.numeric(drug_sig$GeneID)
  
  save(drug_sig, file="GSE69845_drug_MCF7.RData")

  sort(cor(drug_sig[, -1])["Vinblastine",])
  