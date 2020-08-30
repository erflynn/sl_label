
require('tidyverse')
require('glmnet')
require('bestNormalize')

args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
ds <- "sex"
load(file=sprintf("data/07_model_dat/fit_%s_rnaseq_%s.RData", prefix, ds))

my_rows <- unique(rownames(fit$beta))
slStudy <- function(study_id){
  my_dat <- read_csv(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id))  %>% as.data.frame()
  rownames(my_dat) <- my_dat$gene_name
  # select genes and rotate
  return(my_dat[my_rows,] %>% select(-gene_name)) # return labels
}

exp_samp <- read_csv("data/missing_rnaseq_pred.csv") %>% 
  filter(organism==prefix)
list.studies <- unique(exp_samp$study_acc)

f_exist <- sapply(list.studies, function(study_id) 
  file.exists(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id)))
list.studies[!f_exist]

## fill in missing files
sapply(list.studies[!f_exist], function(study_id){
  list_samples <- exp_samp %>% filter(study_acc==study_id) %>% pull(sample_acc)
  chunk_df <- extractStudyChunk(study_id, list_samples) # from 02../01_quant_to_csv.R
  chunk_df %>% write_csv(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id))
})


##



sample_df <- exp_samp %>% 
  filter(present & num_reads >= 100000 & study_acc %in% list.studies)

print("extracting")
all_dat0 <- do.call(cbind, lapply(list.studies, slStudy))
all_dat <- all_dat0[,intersect(colnames(all_dat0), sample_df$sample_acc)] # filter for appropriate ones
all_dat[is.na(all_dat)] <- 0
num_zeros <- apply(all_dat, 1, function(x) sum(x==0))
all_dat2 <- all_dat[num_zeros <= 0.9*ncol(all_dat),]

print("normalizing")
x_test <- apply(all_dat2, 1, function(row) boxcox(row+0.5)$x.t)
x_test2 <- cbind(x_test, t(all_dat[num_zeros > 0.9*ncol(all_dat),]))
print("predicting")
preds <- predict(fit, newx=x_test2[,my_rows], s=my.lambda, type="response")
preds2 <- data.frame("id"=rownames(preds), "pred"=preds[,1])

preds2 %>% 
  write_csv(sprintf("data/%s_rnaseq_sl_m.csv", prefix))
