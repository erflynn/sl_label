require('glmnet')
require('tidyverse')
require('data.table')

options(stringAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

expr_cl <- fread(sprintf("data/microarray/%s/04_sl_input/%s_cell_pos_expr.csv", prefix, prefix ), data.table=FALSE)
expr_no_cl <- fread(sprintf("data/microarray/%s/04_sl_input/%s_cell_neg_expr.csv", prefix, prefix ), data.table=FALSE)

mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample.csv", prefix), data.table=FALSE)

study_no_cl <- mapping %>% filter(sample_acc %in% colnames(expr_no_cl) )
study_cl <- mapping %>% filter(sample_acc %in% colnames(expr_cl))

train_test <- function(df) {
  df2 <- df %>% 
    group_by(study_acc) %>% 
    count() %>% 
    arrange(desc(n)) %>% 
    ungroup() %>% 
    group_by(n) %>% 
    mutate(num_studies=n()) %>% 
    sample_n(num_studies) %>% 
    mutate(ord=1:n()) %>% 
    arrange(desc(n, ord)) %>%
    ungroup() %>%
    mutate(row_num=1:n()) 
    
  test <- df2 %>% 
    filter(row_num %%2==0) %>% 
    select(-n, -num_studies, -ord, -row_num) %>%
    left_join(df, by="study_acc")
  train <- df2 %>% 
    filter(row_num %%2==1) %>% 
    select(-n, -num_studies, -ord, -row_num) %>%
    left_join(df, by="study_acc")
  
  # ensure absolutely no overlap or replicates
  test2 <- test %>% 
    anti_join(train, by="sample_acc") %>%  
    ungroup() %>% select(sample_acc) %>% unique()
  train2 <- train %>% ungroup() %>% select(sample_acc) %>% unique()
  
  return(list("train"=train, "test"=test2))
}

get_train_folds <- function(df, ngroups){
  
  study_counts_df <- df %>% 
    group_by(study_acc) %>% 
    count() %>% 
    arrange(desc(n)) %>%
    ungroup() %>% 
    rename(num_samples=n) %>% 
    group_by(num_samples) %>% 
    mutate(n=n() ) %>%
    sample_n(n) %>% 
    mutate(ord=1:n()) %>% 
    arrange(desc(num_samples, ord))
  
  study_counts <- study_counts_df %>% 
    select(study_acc, num_samples) %>% 
    ungroup() %>%
    filter(!is.na(study_acc))
  
  group_l <- c(1:ngroups, ngroups:1)
  reps <- nrow(study_counts) %/% length(group_l)
  rems <- nrow(study_counts) %% length(group_l)
  opp_group_l <- c(ngroups:1,1:ngroups)
  study_counts$fold <- c(rep(group_l, reps), opp_group_l[1:rems])
  sample_df <- study_counts %>% select(-num_samples) %>% left_join(df, by=c("study_acc"))
  return(sample_df)
}



# divide into training + testing
cl_div <- train_test(study_cl)
no_cl_div <- train_test(study_no_cl)

cl_train <- unique(cl_div$train$sample_acc)
cl_test <- unique(cl_div$test$sample_acc)
no_cl_train <- unique(no_cl_div$train$sample_acc)
no_cl_test <- unique(no_cl_div$test$sample_acc)

train_lab <- c(rep(1, length(cl_train)), rep(0, length(no_cl_train)))
test_lab <- c(rep(1, length(cl_test)), rep(0, length(no_cl_test)))
expr_train <- cbind(expr_cl[,cl_train], expr_no_cl[,no_cl_train])
expr_test <- cbind(expr_cl[,cl_test], expr_no_cl[,no_cl_test])
rownames(expr_train) <- expr_no_cl$rid
rownames(expr_test) <- expr_no_cl$rid

# then CV folds
train_study <- mapping %>% 
  filter(sample_acc %in% colnames(expr_train)) %>%
  mutate(class=ifelse(sample_acc %in% colnames(expr_cl), 1, 0))
train_folds <- get_train_folds(train_study, 10)
train_folds %>% select(-study_acc) %>% group_by(fold, class) %>% count() # check balance
train_folds2 <- train_folds %>% select(sample_acc, fold) %>% unique()
fold_list <- train_folds2$fold
names(fold_list) <- train_folds2$sample_acc
fold_list2 <- fold_list[colnames(expr_train)]
names(fold_list2) <- NULL

# train lasso on these data
x_train <- as.matrix(t(expr_train))
#x_train <- apply(x_train, c(1,2),as.numeric)

x_test <- as.matrix(t(expr_test))
#x_test <- apply(x_test, c(1,2),as.numeric)
y <-factor(train_lab)
x <- x_train



run_fold <- function(my.fold, nfolds){
  train_data <- train_valid %>% filter(!fold %in% my.fold)
  valid_data <- train_valid %>% filter(fold %in% my.fold)
  # redo the training fold
  my.folds <- data.frame(cbind("fold"=unique(train_data$fold), 
                               "new_fold"=1:(nfolds-length(my.fold))))
  train_data2 <- train_data %>% left_join(my.folds)
  
  train_expr_data <- expr_df2[,train_data2$sample_acc]
  train_labels <- train_data2$sex
  train_folds <- train_data2$new_fold
  valid_expr_data <- expr_df2[,valid_data$sample_acc]
  valid_labels <- valid_data$sex
  
  # transform and run 
  x_train <- t(train_expr_data)
  x_valid <- t(valid_expr_data)
  
  x_train2 <- apply(x_train, 2, function(col) boxcox(col+0.5)$x.t)
  x_valid2 <- apply(x_valid, 2, function(col) boxcox(col+0.5)$x.t)
  
  test_params <- function(my.alpha){
    cvfit = cv.glmnet(x_train2, train_labels, 
                      family="binomial",
                      foldid=train_folds,
                      nfolds=5,
                      alpha=my.alpha, 
                      standardize=FALSE)
    preds_train <- predict(cvfit, newx=x_train2, s="lambda.1se", type="response")
    preds_class_train <- sapply(predict(cvfit, newx=x_train2, s="lambda.1se", type="class"), 
                                as.numeric)
    train_acc <- sum(preds_class_train==train_labels)/length(train_labels)
    print(train_acc)
    preds_class_valid <- sapply(predict(cvfit, newx=x_valid2, s="lambda.1se", type="class"), as.numeric)
    valid_acc <- sum(preds_class_valid==valid_labels)/length(valid_labels) 
    print(valid_acc)
    return(list("t"=train_acc, 
                "v"=valid_acc, 
                "lambda"=cvfit$lambda.1se, 
                "alpha"=my.alpha))
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


cv.5=cv.glmnet(x,y,foldid=fold_list2, nfolds=10, alpha=.5, 
               family="binomial", 
               standardize=FALSE)
preds_class_train <- sapply(predict(cvfit, newx=x_train, s="lambda.1se", 
                                    type="class"), as.numeric)
train_acc <- sum(preds_class_train==train_lab)/length(train_lab)
print(train_acc)

preds_class_valid <- sapply(predict(cvfit, newx=x_test, s="lambda.1se", type="class"), as.numeric)
valid_acc <- sum(preds_class_valid==test_lab)/length(test_lab) 
print(valid_acc)


cv1=cv.glmnet(x,y,foldid=fold_list2,alpha=1, family="binomial", 
              standardize=FALSE)
cv0=cv.glmnet(x,y,foldid=fold_list2,alpha=0, family="binomial", standardize=FALSE)

#cvfit = cv.glmnet(x_train, train_cl_lab, family="binomial", 
#                 nlambda=200, alpha=0.5, standardize=FALSE, type.measure = "class")
preds_train <- predict(cvfit, newx=x_train, s="lambda.1se", type="response")
preds_class_train <- sapply(predict(cvfit, newx=x_train, s="lambda.1se", type="class"), as.numeric)
sum(preds_class_train==train_cl_lab)/length(train_cl_lab) # 0.986
preds_test <- predict(cvfit, newx=x_test, s="lambda.1se", type="response")
preds_class_test <- sapply(predict(cvfit, newx=x_test, s="lambda.1se", type="class"), as.numeric)
sum(preds_class_test==test_cl_lab)/length(test_cl_lab) # 0.944
mat_coef <- coef(cvfit, lambda="lambda.1se") %>% as.matrix()
save(cvfit, file="data/cvfit_human_cl.RData")
nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef))
load(sprintf("../sex_labeling/geo_pipeline/gpl_ref/%s/%s_gene_map.RData", prefix, prefix))
coef_df2 <- coef_df %>% 
  left_join(ref_dat %>% select("ensembl_gene_id", "hgnc_symbol", "chromosome_name") %>% unique(), 
            by=c("gene"="ensembl_gene_id"))
coef_df2 %>% arrange(coef) %>% write_csv("data/human_cl_coef.csv")
