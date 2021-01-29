require('data.table')
require('tidyverse')
require('groupdata2')
require('glmnet')
require('bestNormalize')
#install.packages("~/Downloads/glmnet_3.0-2.tgz", repos=NULL)
# **NEED GLMNET3.0**



source("code/04_models/cv_utils.R")

options(stringsAsFactors = FALSE)
NFOLDS <- 6
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
data_type <- args[2]
ds <- args[3]
#idx <- as.numeric(args[4])

if (data_type=="microarray"){
  standardizeFlag = TRUE
  
} else {
  standardizeFlag = FALSE
}
# load processed data
pre_processed_dat <- sprintf("data/07_model_dat/%s_%s_%s_train_mat.RData", 
                               prefix, data_type, ds)
  
load(pre_processed_dat) # --> X_train, Y_train, X_test, Y_test


# set up the data
fold_file <- sprintf("data/07_model_dat/%s_%s_fold_breakdown.RData", prefix, data_type)

if (!file.exists(fold_file)){
  # load mapping data
  if (data_type == "microarray"){
    
    mapping <- fread(sprintf("data/01_metadata/%s_compendia_dedup_mapping.csv", 
                             prefix), data.table = FALSE) %>%
      rename(study_acc=study_str)
  } else {
    MIN.READS <- 100000
    # RNA-seq data does not have duplicates... interesting
    mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", 
                             prefix), data.table=FALSE)  %>%
      filter(present & num_reads >= MIN.READS)
  }
  # load processed data #

  pos_samples <- rownames(X_train)[Y_train==1]#, rownames(X_test)[Y_test==1])
  neg_samples <- rownames(X_train)[Y_train==0]#, rownames(X_train)[Y_train==0])
  
  pos_map <- mapping %>% filter(sample_acc %in% pos_samples)
  neg_map <- mapping %>% filter(sample_acc %in% neg_samples)
  
  # deal with class imbalance
  set.seed(304)
  num_pos <- nrow(pos_map)
  num_neg <- nrow(neg_map)
  
  if ((num_pos > (1.2*num_neg)) | (num_neg > (1.2*num_pos))){
    num_ea <- min(num_neg, num_pos)
    pos_map <- pos_map %>% sample_n(num_ea)
    neg_map <- neg_map %>% sample_n(num_ea)
    
  }
  
  my_dat <- bind_rows(pos_map %>% mutate(class=1), 
                      neg_map %>% mutate(class=0)) %>%
    mutate(class=factor(class), study_acc=factor(study_acc))
  
  # training + testing division
  #if (nrow(my_dat) < 1000){
  #  tfolds <- 4
  #  NFOLDS <- 4
  #} else {
  #  tfolds <- 5
  #}
  
  train_test_split <- partition_group_data(my_dat, grp_col="study_acc", 
                                           nfolds=5)
  # train_test_split %>% group_by(partition, class) %>% count()
  # train_part <- train_test_split %>% 
  #   filter(partition!=tfolds) %>% 
  #   arrange(class, sample_acc)
  # test_part <- train_test_split %>% 
  #   filter(partition==tfolds) %>% 
  #   arrange(class, sample_acc)
  # 
  # 
  # divide into folds
  fold_list <- partition_group_data(my_dat, grp_col="study_acc",
                                    nfolds=NFOLDS) %>% 
    arrange(class, sample_acc)
  samp_to_fold <- fold_list$partition
  names(samp_to_fold) <- fold_list$sample_acc
  save(fold_list, file=fold_file)
}
load(fold_file)





# run a CV fold with a particular alpha
test_params <- function(my.alpha, X_train2, Y_train2, X_valid2, Y_valid2, train_folds){
  cvfit = cv.glmnet(X_train2, Y_train2, 
                    family="binomial",
                    foldid=train_folds,
                    alpha=my.alpha, 
                    standardize=standardizeFlag,
                    trace=FALSE)
  
  train_assess <- assess.glmnet(cvfit, newx=X_train2, newy=Y_train2)
  valid_pred <- predict(cvfit, newx=X_valid2, type="response")
  valid_assess <- assess.glmnet(cvfit, newx=X_valid2, newy=Y_valid2)
  my.lambda <- cvfit$lambda.1se
  
  ta <- unlist(train_assess) 
  va <- unlist(valid_assess) 
  train_valid_assess <- bind_rows(ta, va)
  train_valid_assess$grp <- c("train", "valid")
  train_valid_assess$alpha <- my.alpha
#  train_valid_assess$ngenes <- ngenes
  train_valid_assess$lambda <- my.lambda
  
  # turn this into a dataframe
  # fold, my.alpha, adaptive, lambda, train_assess, valid_assess
  
  # inclue valid_pred, Y_valid
  df <- data.frame(cbind("sample"=rownames(valid_pred), 
                         "yhat"=valid_pred[,1]))
  df$y <- Y_valid2
  rownames(df) <- NULL
  
  return(list("tv"=train_valid_assess, "ydf"=df))
}

# run a full CV fold
run_fold <- function(my.fold, nfolds){
  train_data <- fold_list %>% filter(!partition %in% my.fold) %>% arrange(partition)
  valid_data <- fold_list %>% filter(partition %in% my.fold)
  # redo the training fold
  my.folds <- data.frame(cbind("partition"=unique(train_data$partition), 
                               "new_fold"=1:(nfolds-length(my.fold))))

  train_data2 <- train_data %>% left_join(my.folds)
  
  X_train1 <- X_train[train_data2$sample_acc,]
  Y_train2 <- train_data2$class
  train_folds <- train_data2$new_fold
  X_valid1 <- X_train[valid_data$sample_acc,]
  Y_valid2 <- valid_data$class


  res1 <- lapply(seq(0,1,0.1), function(my.alpha) 
    test_params(my.alpha, X_train1, Y_train2, X_valid1, Y_valid2, train_folds) )
  return(res1)
}

# run everything
full_res <- lapply(1:(NFOLDS), function(idx) {
  print(idx)
  res <- run_fold(idx, NFOLDS)
  save(res, file=sprintf("data/06_fold_dat/fold_%s_%s_%s_%s.RData", 
                         prefix, data_type, ds, idx))
  
  return(res)
})
