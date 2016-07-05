setwd("~/Documents/stanford/tumor_cell_line/RGES_manuscript/release/data")

cancers = c("BRCA", "LIHC", "COAD")

#get disease signatures
source("../code/disease_sig.R")

#compute RGES of compounds with drug efficacy data. Either landmark genes or whole-genome can be used.
#support other methods to measure reverse potency (spearman, pearson, cosine)
source("../code/RGES_computation_cmpd_selected.R")

#compute RGES of all compounds using the landmark genes
source("../code/RGES_computation_cmpd_all.R")

#lincs statistics
source("../code/lincs_stats.R")

#confounding factors
source("../code/RGES_confounding.R")

#dose vs RGES
source("../code/RGES_dose_cor.R")

#RGES vs connectivity cmap
source("../code/RGES_cmap_score_diff.R")

#RGES in single cell lines
source("../code/RGES_IC50_in_single_cell_lines.R")

#explain outliers
source("../code/RGES_vinblastine_external.R")

#sRGES vs IC50 
source("../code/sRGES_IC50.R")

#assessing performance using different methods
source("../code/performance.R")

#use sRGES to prioritize compounds
source("../code/sRGES_all_cmpds.R")

#sRGES vs IC50 in LIHC after adding validation data
source("../code/sRGES_IC50_LIHC_after_validation.R")

#apply sRGES method to the connectivity score
#comparing with the summarization method proposed in LINCS
source("../code/sRGES_connectivity_score_comparison.R")

#identify reversed genes
source("../code/compute_reversed_genes.R")
