require('data.table')
require('tidyverse')
require('groupdata2')
require('glmnet')
#install.packages("~/Downloads/glmnet_3.0-2.tgz", repos=NULL)
# **NEED GLMNET3.0**

source("code/03_models/cv_utils.R")

options(stringsAsFactors = FALSE)
NFOLDS <- 8
prefix <- "human"
data_type <- "microarray"
ds <- "cell"
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
data_type <- args[2]
ds <- args[3]
i <- as.numeric(args[4])

# set up the data
pos <- fread(sprintf("data/%s/%s/04_sl_input/%s_%s_pos_expr.csv", data_type, prefix, prefix, ds ), data.table=FALSE)
neg <- fread(sprintf("data/%s/%s/04_sl_input/%s_%s_neg_expr.csv", data_type, prefix, prefix, ds ), data.table=FALSE)

if (data_type == "microarray"){
  #mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample.csv", prefix), data.table=FALSE)  
  mapping <- fread(sprintf("data/01_metadata/%s_compendia_dedup_mapping.csv", prefix, data_type), data.table = FALSE) %>%
    rename(study_acc=study_str)
} else {
  MIN.READS <- 100000
  # RNA-seq data does not have duplicates... interesting
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix), data.table=FALSE)  %>%
    filter(present & num_reads >= MIN.READS)
}



pos_map <- mapping %>% filter(sample_acc %in% colnames(pos))
neg_map <- mapping %>% filter(sample_acc %in% colnames(neg))
my_dat <- bind_rows(pos_map %>% mutate(class=1), 
                    neg_map %>% mutate(class=0)) %>%
  mutate(class=factor(class), study_acc=factor(study_acc))

# training + testing division
set.seed(304)
train_test_split <- partition_group_data(my_dat, grp_col="study_acc", nfolds=5)
train_test_split %>% group_by(partition, class) %>% count()
train_part <- train_test_split %>% filter(partition!=5) %>% arrange(class, sample_acc)
test_part <- train_test_split %>% filter(partition==5) %>% arrange(class, sample_acc)

Y_train <- train_part$class
Y_test <- test_part$class

# divide into folds
fold_list <- partition_group_data(train_part %>% select(-partition), nfolds=NFOLDS) %>% arrange(class, sample_acc)
samp_to_fold <- fold_list$partition
names(samp_to_fold) <- fold_list$sample_acc
  
# grab relevant dfs
X_train <- t(cbind(neg[,train_part$sample_acc[Y_train==0]],pos[,train_part$sample_acc[Y_train==1]]))
X_test <- t(cbind(neg[,test_part$sample_acc[Y_test==0]],pos[,test_part$sample_acc[Y_test==1]]))
colnames(X_train) <- neg$rid
colnames(X_test) <- neg$rid
rownames(X_train) <- c(train_part$sample_acc)
rownames(X_test) <- c(test_part$sample_acc)


# --- preprocess as needed (e.g. boxcox for RNASeq) --- #
if (data_type=="rnaseq"){
  
  X_train <- apply(X_train, 2, function(col) boxcox(col+0.5)$x.t) 
  X_test <- apply(X_test, 2, function(col) boxcox(col+0.5)$x.t)
  
}

# CV folds
rand_genes <- sample(1:ncol(X_test), 2000)


run_fold <- function(my.fold, nfolds=NFOLDS){
  train_data <- fold_list %>% filter(!partition %in% my.fold)
  valid_data <- fold_list %>% filter(partition %in% my.fold)
  # redo the training fold
  my.folds <- data.frame(cbind("partition"=unique(train_data$partition), "new_fold"=1:(nfolds-length(my.fold))))
  train_data2 <- train_data %>% left_join(my.folds)
  
  X_train2 <- X_train[train_data2$sample_acc,]
  Y_train2 <- train_data2$class
  train_folds <- train_data2$new_fold
  X_valid <- X_train[valid_data$sample_acc,]
  Y_valid <- valid_data$class
  
  if (data_type=="microarray"){
    standardizeFlag = TRUE
  } else {
    standardizeFlag = FALSE
  }
  
  test_params <- function(my.alpha, adaptive, relax){
    # add adaptive lasso
    if (adaptive){
      ridge_cv = cv.glmnet(X_train2[,rand_genes], Y_train2, 
                        family="binomial",
                        foldid=train_folds,
                        alpha=0, 
                        standardize=standardizeFlag,
                        trace=TRUE)
      best_ridge_coef <- coef(ridge_cv, s = ridge_cv$lambda.min)
      best_ridge_coef <- as.numeric(best_ridge_coef)[-1]
      
      cvfit = cv.glmnet(X_train2[,rand_genes], Y_train2, 
                        family="binomial",
                        foldid=train_folds,
                        alpha=my.alpha, 
                        standardize=standardizeFlag,
                        trace=TRUE,
                        penalty.factor = 1 / abs(best_ridge_coef))
    } else {
      cvfit = cv.glmnet(X_train2[,rand_genes], Y_train2, 
                        family="binomial",
                        foldid=train_folds,
                        alpha=my.alpha, 
                        standardize=standardizeFlag,
                        trace=TRUE)
    }
    
    train_assess <- assess.glmnet(cvfit, newx=X_train2[,rand_genes], newy=Y_train2)
    valid_pred <- predict(cvfit, newx=X_valid[,rand_genes], type="response")
    valid_assess <- assess.glmnet(cvfit, newx=X_valid[,rand_genes], newy=Y_valid)
    my.lambda <- cvfit$lambda.1se
    
    ta <- unlist(train_assess) 
    va <- unlist(valid_assess) 
    train_valid_assess <- bind_rows(ta, va)
    train_valid_assess$grp <- c("train", "valid")
    train_valid_assess$alpha <- my.alpha
    train_valid_assess$adaptive <- adaptive
    train_valid_assess$lambda <- my.lambda
    
    # turn this into a dataframe
    # fold, my.alpha, adaptive, lambda, train_assess, valid_assess
    
    # inclue valid_pred, Y_valid
    df <- data.frame(cbind("sample"=rownames(valid_pred), "yhat"=valid_pred[,1]))
    df$y <- Y_valid
    rownames(df) <- NULL
    
    return(list("tv"=train_valid_assess, "ydf"=df))
  }
  
  param_res <- data.frame()
  for (my.alpha in seq(0,1, 0.1)){
    alpha_res <- test_params(my.alpha)
    if (my.alpha==0){
      param_res <- alpha_res
    } else {
      param_res <- rbind(param_res, alpha_res)
    }
  }
  return(param_res)
}

res <- run_fold(folds[i])
write_csv(res, file=sprintf("data/%s/%s/04_sl_input/fold_%s_%s.csv", data_type, prefix, ds, i))

# fold_res <- data.frame()
# folds <- unique(train_valid$fold)
# for (i in 1:length(folds)){
#   print(i)
#   res <- run_fold(folds[i])
#   if (i==1){
#     fold_res <- res
#   }
#   else {
#     fold_res <- rbind(fold_res, res)
#   }
# } 



# for each fold in folds
#  for each param in params:
#    cv.glmnet() on all but fold using params (alpha)
#    predict on fold
#    return(pred_accuracy)