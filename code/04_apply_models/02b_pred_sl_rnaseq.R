require('tidyverse')
require('glmnet')
require('bestNormalize')

args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
idx <- as.numeric(args[2])
ds <- "sex"
extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}

# prediction function, list of genes
load("data/07_model_dat/human_rnaseq_boxcox.RData") # --> list_genes, box_cox_obj
load(file=sprintf("data/07_model_dat/fit_%s_rnaseq_%s.RData", prefix, ds))
my_rows <- unique(rownames(fit$beta))
stopifnot(my_rows==list_genes)
# slStudy <- function(study_id){
#   my_dat <- read_csv(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id))  %>% as.data.frame()
#   rownames(my_dat) <- my_dat$gene_name
#   my_dat$gene_name <- NULL 
#   # select genes and rotate
#   return(my_dat[list_genes,]) # return labels
# }
slStudy <- function(study_id){
  my_dat <- read_csv(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id))  %>% as.data.frame()
  rownames(my_dat) <- my_dat$gene_name
  # select genes and rotate
  return(my_dat[my_rows,] %>% select(-gene_name)) # return labels
}


exp_samp <- read_csv(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix))
list.studies <- unique(exp_samp$study_acc)
f_exist <- sapply(list.studies, function(study_id) 
  file.exists(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id)))
sample_df <- exp_samp %>% 
  filter(present & study_acc %in% list.studies[f_exist])

study_chunk <- extractChunk(list.studies[f_exist], idx, 50)
print("extracting")
all_dat0 <- do.call(cbind, lapply(study_chunk, slStudy))
all_dat <- t(all_dat0[,intersect(colnames(all_dat0),sample_df$sample_acc)]) # filter for appropriate ones

print("normalizing")
all_dat2 <- do.call(cbind, lapply(1:length(list_genes), 
                                 function(i) predict(box_cox_obj[[i]], 
                                                     newdata=(all_dat[,i]+0.5) )
))
#x_test <- apply(all_dat2, 1, function(row) boxcox(row+0.5)$x.t)
#x_test2 <- cbind(x_test, t(all_dat[num_zeros > 0.9*ncol(all_dat),]))[,my_rows]
print("predicting")
preds <- predict(fit, newx=all_dat2, s=my.lambda, type="response")
preds2 <- data.frame("id"=rownames(all_dat), "pred"=preds[,1])
  
preds2 %>% write_csv(sprintf("data/08_model_out/%s_rnaseq_%s_%s_sl.csv", prefix, ds, idx))
