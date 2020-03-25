require('tidyverse')
require('glmnet')
require('bestNormalize')

args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
idx <- as.numeric(args[2])

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}

# prediction function, list of genes
load(file=sprintf("data/rnaseq/%s/04_model_out/rnaseq_sl_fit.RData", prefix)) # my.lambda, fit
my_rows <- unique(rownames(fit$beta))
slStudy <- function(study_id){
  my_dat <- read_csv(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id))  %>% as.data.frame()
  rownames(my_dat) <- my_dat$gene_name
  my_dat$gene_name <- NULL 
  # select genes and rotate
  return(my_dat[my_rows,]) # return labels
}

exp_samp <- read_csv(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix))
list.studies <- unique(exp_samp$study_acc)
sample_df <- exp_samp %>% 
  filter(present & num_reads >= 100000 & study_acc %in% list.studies)

study_chunk <- extractChunk(list.studies, idx, 50)
print("extracting")
all_dat0 <- do.call(cbind, lapply(study_chunk, slStudy))
all_dat <- all_dat0[,sample_df$sample_acc] # filter for appropriate ones
all_dat[is.na(all_dat)] <- 0
num_zeros <- apply(all_dat, 1, function(x) sum(x==0))
all_dat2 <- all_dat[num_zeros <= 0.9*ncol(all_dat),]

print("normalizing")
x_test <- apply(all_dat2, 1, function(row) boxcox(row+0.5)$x.t)
x_test2 <- cbind(x_test, t(all_dat[num_zeros > 0.9*ncol(all_dat),]))[,my_rows]
print("predicting")
preds <- predict(fit, newx=x_test2, s=my.lambda, type="response")
preds2 <- data.frame("id"=rownames(preds), "pred"=preds[,1])
  
preds2 %>% write_csv(sprintf("data/rnaseq/%s/04_model_out/%s_rnaseq_sl_%s.csv", prefix, prefix, idx))
