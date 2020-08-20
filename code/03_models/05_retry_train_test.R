
# --- update the RNA-seq 

my_files =list.files(pattern="fold_mouse_rnaseq_sex*")
res_list <- do.call(rbind, lapply(my_files, function(my_f) {
  load(my_f); 
  res_df <- do.call(rbind, lapply(res, function(x) 
    do.call(rbind, lapply(x, function(y) y$tv))))
  return(res_df)
  }))
res_list %>% filter(grp=="valid") %>% arrange(class.1)


res_list %>% filter(grp=="valid") %>% 
  group_by(alpha, ngenes) %>% 
  summarize(across(c(deviance.1:mae.1, lambda), median)) %>%
  arrange(class.1)

# this actually looks good! esp for RNA-seq data. YAYYYY
# --- for the cell line data not so much sadly --- #

# ---------------------------- #

# play with training and testing data...
require('glmnet')

load("data/07_model_dat/human_rnaseq_sex_train_mat.RData")
# --> X_test, X_train
# --> Y_test, Y_train
require('limma')

fit <- lmFit(t(X_train), as.numeric(Y_train))
fit <- eBayes(fit)
tt <-topTable(fit, number=ncol(X_train))
tt$transcript <- rownames(tt)
tt %>% head(NUM_GENES)
# standardize == FALSE?