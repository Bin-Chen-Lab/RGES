setwd("~/Documents/stanford/tumor_cell_line/RGES_manuscript/release/data")

code_dir = "../git_code/RGES/"
cancers = c("BRCA", "LIHC", "COAD", "ER")

#core functions
source(paste(code_dir, "core_functions.R", sep=""))

#get disease signatures
source(paste(code_dir, "disease_sig.R", sep=""))

#compute RGES of compounds with drug efficacy data. Either landmark genes or whole-genome can be used.
#support other methods to measure reverse potency (spearman, pearson, cosine)
source(paste(code_dir, "RGES_computation_cmpd_selected.R", sep=""))

#compute RGES of all compounds using the landmark genes
source(paste(code_dir, "RGES_computation_cmpd_all.R", sep=""))

#LINCS statistics
source(paste(code_dir, "lincs_stats.R", sep=""))

#confounding factors. RGES depends on cell line, treatment duration, treatment concentration
source(paste(code_dir, "RGES_confounding.R", sep=""))

#relationship between treatment duration vs RGES
source(paste(code_dir, "RGES_dose_cor.R", sep=""))

#Differences between RGES and connectivity cmap score. 
source(paste(code_dir, "RGES_cmap_score_diff.R", sep=""))

#correlation between RGES and IC50 within single cell lines
source(paste(code_dir, "RGES_IC50_in_single_cell_lines.R", sep=""))

#explain outliers (which do not follow the linear correlation)
source(paste(code_dir, "RGES_vinblastine_external_prepare.R", sep=""))
source(paste(code_dir, "RGES_vinblastine_external_analysis.R", sep=""))

#sRGES vs IC50 
source(paste(code_dir, "sRGES_IC50.R", sep=""))

#assessing performance using different methods
source(paste(code_dir, "performance.R", sep=""))

#use sRGES to prioritize compounds
source(paste(code_dir, "sRGES_all_cmpds.R", sep=""))

#sRGES vs IC50 in LIHC after adding validation data
source(paste(code_dir, "sRGES_IC50_LIHC_after_validation.R", sep=""))

#apply sRGES method to the connectivity score
#comparing with the summarization method proposed in LINCS
source(paste(code_dir, "sRGES_connectivity_score_comparison.R", sep=""))

#identify reversed genes
source(paste(code_dir, "compute_reversed_genes_v2=.R", sep=""))
