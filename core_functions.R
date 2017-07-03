##############
cmap_score_new <- function(sig_up, sig_down, drug_signature) {
  #the old function does not support the input list with either all up genes or all down genes, this new function attempts to addess this.
  #we also modify the original CMap approach: whenever the sign of ks_up/ks_down, we substract the two scores such that the final scores would not enrich at 0.
  
  num_genes <- nrow(drug_signature)
  ks_up <- 0
  ks_down <- 0
  connectivity_score <- 0
  
  # I think we are re-ranking because the GeneID mapping changed the original rank range
  drug_signature[,"rank"] <- rank(drug_signature[,"rank"])
  
  # Merge the drug signature with the disease signature by GeneID. This becomes the V(j) from the algorithm description
  up_tags_rank <- merge(drug_signature, sig_up, by.x = "ids", by.y = 1)
  down_tags_rank <- merge(drug_signature, sig_down, by.x = "ids", by.y = 1)
  
  up_tags_position <- sort(up_tags_rank$rank)
  down_tags_position <- sort(down_tags_rank$rank)
  
  num_tags_up <- length(up_tags_position)
  num_tags_down <- length(down_tags_position)
  
  # 
  if(num_tags_up > 1) {
    a_up <- 0
    b_up <- 0
    
    a_up <- max(sapply(1:num_tags_up,function(j) {
      j/num_tags_up - up_tags_position[j]/num_genes
    }))
    b_up <- max(sapply(1:num_tags_up,function(j) {
      up_tags_position[j]/num_genes - (j-1)/num_tags_up
    }))
    
    if(a_up > b_up) {
      ks_up <- a_up
    } else {
      ks_up <- -b_up
    }
  }else{
    ks_up <- 0
  }
  
  if (num_tags_down > 1){
    
    a_down <- 0
    b_down <- 0
    
    a_down <- max(sapply(1:num_tags_down,function(j) {
      j/num_tags_down - down_tags_position[j]/num_genes
    }))
    b_down <- max(sapply(1:num_tags_down,function(j) {
      down_tags_position[j]/num_genes - (j-1)/num_tags_down
    }))
    
    if(a_down > b_down) {
      ks_down <- a_down
    } else {
      ks_down <- -b_down
    }
  }else{
    ks_down <- 0
  }
  
  if (ks_up == 0 & ks_down != 0){ #only down gene inputed
    connectivity_score <- -ks_down 
  }else if (ks_up !=0 & ks_down == 0){ #only up gene inputed
    connectivity_score <- ks_up
  }else if (sum(sign(c(ks_down,ks_up))) == 0) {
    connectivity_score <- ks_up - ks_down # different signs
  }else{
    connectivity_score <- ks_up - ks_down
  }
  
  return(connectivity_score)
}

get.gene.list <- function(con){
  lincs_gene <- dbReadTable(con, "probe_id_info")
  
  return(lincs_gene$gene_id)
}

get.instance.sig <- function(id, con,  landmark=F){
  #
  
  sig_file <- paste("~/Documents/stanford/lincs/data/lincs/", id, ".txt", sep="")
  sig_value <- NULL
  
  if (file.exists(sig_file)){
    sig_value <- scan(sig_file )
    
    if (landmark){
      sig_value <- sig_value[1:978]
    }
    
  }else{
    query <- paste("select * from proj_lincs.sig_values where id = ", id, sep="")
    rs <- dbSendQuery(con, query)
    result <- fetch(rs, n = -1)[1,]
    value <- as.double(unlist(strsplit(result[1,2], ",")))
    if (landmark == T){
      instance.sig <- data.frame(sig_id =  id, probe_id = seq(1, 978), value[1:978])      
    }else{
      instance.sig <- data.frame(sig_id = id, probe_id = seq(1, length(value)), value = value)      
    }
    dbClearResult(rs)  
    sig_value <- instance.sig$value
  }
  return (sig_value)
}


#find best alpha and beta
find_alpha_beta <- function(){
  alphas <- seq(-1, 1, 0.1)
  betas <- seq(-1, 1, 0.1)
  all_values <- data.frame()
  for (alpha in alphas){
    for (beta in betas){
      lincs_drug_prediction_subset <- subset(lincs_drug_prediction, cell_id %in% c(cell_line_selected)) #HT29 MCF7
      lincs_drug_prediction_subset$RGES <- sapply(1:nrow(lincs_drug_prediction_subset), function(id){
        getsRGES(lincs_drug_prediction_subset[id,"RGES"], lincs_drug_prediction_subset[id, "pert_dose"], lincs_drug_prediction_subset[id, "pert_time"], alpha, beta)
      })
      lincs_drug_prediction_subset <- aggregate(RGES ~ pert_iname, lincs_drug_prediction_subset, mean)
      
      activity_RGES <- merge(lincs_drug_prediction_subset, lincs_drug_activity_subset, by="pert_iname")
      
      activity_RGES_summarized <- activity_RGES #aggregate(cbind(RGES, standard_value) ~ pert_iname, activity_RGES,  min)
      
      cor <- cor(activity_RGES_summarized$RGES, log(activity_RGES_summarized$standard_value, 10), method="spearman")
      all_values <- rbind(all_values, data.frame(cor, alpha, beta))
    }
  }
  return(all_values)
}

getsRGES1 <- function(RGES, cor, pert_dose, pert_time){
  sRGES <- RGES
  if (pert_time == 24){
    sRGES <- RGES + predict(lm_dose_24, data.frame(dose=round(log(pert_dose, 10), 1)))
  }
  if (pert_time == 6){
    sRGES <- RGES + predict(lm_dose_6, data.frame(dose=round(log(pert_dose, 10), 1)))
  }
  return (sRGES * cor )
}

getsRGES2 <- function(RGES, cor, pert_dose, pert_time){
  sRGES <- RGES 
  
  #older version
  if (pert_time < 24){
    sRGES <- sRGES - 0.1
  }
  
  if (pert_dose < 10){
    sRGES <- sRGES - 0.2
  }
  return(sRGES * cor)
}

getsRGES3 <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
  sRGES <- RGES
  pert_time <- ifelse(pert_time < 24, "short", "long")
  pert_dose <- ifelse(pert_dose < 10, "low", "high")
  if (pert_time == "short" & pert_dose == "low"){
    sRGES <- sRGES + diff[4]
  }
  if (pert_dose ==  "low" & pert_time == "long"){
    sRGES <- sRGES + diff[2]
  }
  if (pert_dose ==  "high" & pert_time == "short"){
    sRGES <- sRGES + diff[1]
  }
  return(sRGES ) #* cor/max_cor
}

##################
####
getsRGES <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
  sRGES <- RGES
  pert_time <- ifelse(pert_time < 24, "short", "long")
  pert_dose <- ifelse(pert_dose < 10, "low", "high")
  if (pert_time == "short" & pert_dose == "low"){
    sRGES <- sRGES + diff[4]
  }
  if (pert_dose ==  "low" & pert_time == "long"){
    sRGES <- sRGES + diff[2]
  }
  if (pert_dose ==  "high" & pert_time == "short"){
    sRGES <- sRGES + diff[1]
  }
  return(sRGES * cor/max_cor) #
}