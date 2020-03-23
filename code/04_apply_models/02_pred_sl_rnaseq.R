require('tidyverse')
require('glmnet')
args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
idx <- as.numeric(args[2])

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}

# prediction function, list of genes
load(file=sprintf("data/%s/04_sl_input/rnaseq_fun.RData", prefix)) # pred_fun
xy_genes <- read_csv(sprintf("data/%s/04_sl_input/xy_genes_rnaseq.csv", prefix)) 

slStudy <- function(study_id){
  my_dat <- read_csv(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id))  %>% as.data.frame()
  rownames(my_dat) <- "gene_name"
  my_dat$gene_name <- NULL 
  # select genes and rotate
  my_dat2 <- my_dat[xy_genes,] # select the genes of interest

  preds <- pred_fun(t(my_dat2), type="response") #predict(fit, newx=t(my_dat), s="lambda.1se", type="response")
  preds2 <- data.frame("id"=rownames(preds), "pred"=preds[,1])
  # label
  
  return(preds2) # return labels
}

preds2 <- 

preds2 %>% write_csv(sprintf("data/%s/05_sl_output/%s_rnaseq_sl_%s.csv", prefix, prefix, idx))
