require('tidyverse')
require('glmnet')
require('bestNormalize')
require('limma')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]



## ------- SETUP -------- ##

load(sprintf("data/rnaseq/%s/03_model_in/expr_sl_in.RData", prefix)) # --> all_df

sl <- read_csv(sprintf("data/01_metadata/%s_rnaseq_sex_lab.csv",prefix)) %>% 
  as.data.frame()
rownames(sl) <- sl$sample_acc



# do something different for rat
if (prefix == "rat") {
  samples <- intersect(colnames(all_df), sl$sample_acc)
  expr_df <- all_df[,samples]
  meta_df <- sl[samples,]
} else {
  # filter so no cell line samples (?)
  metadata <- read.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", prefix))
  metadata_sm <- metadata %>% filter(acc %in% colnames(all_df))
  no_cl <- metadata_sm %>% 
    filter(cl_line %in% c("", "--", "n/a", "not applicable", "na"))
  
  samples <- intersect(no_cl$acc, sl$sample_acc)
  expr_df <- all_df[,samples]
  meta_df <- sl[samples,] 
}



# setup genes
xy_genes <- read_csv(sprintf("data/rnaseq/%s/03_model_in/xy_genes_rnaseq.csv", prefix))
rownames(expr_df) <- all_df$gene_name
expr_df2 <- expr_df[xy_genes$transcript,]

# remove mostly empty rows
expr_df2[is.na(expr_df2)] <- 0
num_zeros <- apply(expr_df2, 1, function(x) sum(x==0)) # remove > 70% zeros
expr_df2 <- expr_df2[(num_zeros/ncol(expr_df2) <= 0.7),]


train_fold <- function(fold_df, expr_df4, fold_i){
  train_dat <- fold_df %>% filter(fold!=fold_i)
  test_dat <- fold_df %>% filter(fold==fold_i)
  
  design <- train_dat$sex; exp <- expr_df4[train_dat$sample_acc,]
  fit <- lmFit(t(exp), design)
  fit <- eBayes(fit)
  tt <-topTable(fit, number=ncol(exp))
  tt$transcript <- rownames(tt)
  top_f <- tt %>% arrange(logFC) %>% head()
  top_m <- tt %>% arrange(desc(logFC)) %>% head()
  
  list_transcripts <- c(top_f$transcript, top_m$transcript)
  
  # BUILD A MODEL...
  fit = glmnet(expr_df4[train_dat$sample_acc,list_transcripts], train_dat$sex, 
               family="binomial",
               alpha=1, 
               standardize=FALSE)
  preds_class_train <- sapply(predict(fit, newx=expr_df4[train_dat$sample_acc,list_transcripts],
                                      s=0.00208, type="class"),  as.numeric)
  train_acc <- sum(preds_class_train==train_dat$sex)/length(train_dat$sex)
  print(train_acc)
  table(data.frame(cbind(preds_class_train, train_dat$sex)))
  
  #preds_valid <- predict(cvfit, newx=x_valid, s="lambda.1se", type="response")
  preds_class_valid <- sapply(predict(fit, newx=expr_df4[test_dat$sample_acc,list_transcripts], 
                                      s=0.00208, type="class"), as.numeric)
  test_acc <- sum(preds_class_valid==test_dat$sex)/length(test_dat$sex) 
  print(test_acc)
  table(data.frame(cbind(preds_class_valid, test_dat$sex)))
  
  return(list("fold"=fold_i, "pl"=preds_class_valid,"tl"=test_dat$sex))
}

if (prefix == "rat"){
  expr_df2$gene_name <- NULL
  expr_df3 <- apply(expr_df2, 1, function(row) boxcox(row+0.5)$x.t)
  missing <- apply(expr_df3, 2, function(x) sum(is.nan(x)))
  expr_df4 <- expr_df3[,which(missing == 0)]
  study_counts_sex <- meta_df %>% filter(sample_acc %in% colnames(expr_df2)) %>% group_by(study_acc) %>%
    summarize(num_samples=n(), num_f=sum(sex=="female"), num_m=sum(sex=="male")) %>% arrange(desc(num_f), desc(num_m))
  study_counts_sex$fold <- c(rep(c(1, 2, 3, 3, 2, 1), 2), c(3,2))
  study_counts_sex %>% ungroup() %>% group_by(fold) %>%
    select(-study_acc) %>%
    summarize_all(sum)
  meta_df2 <- meta_df %>% mutate(sex=ifelse(sex=="female", 0, 1)) 
  fold_df <- study_counts_sex %>% select(study_acc, fold) %>% left_join(meta_df2) 
  
  # SELECT PREDICTORS
  fold_res <- lapply(1:3, function(i) train_fold(fold_df, expr_df4, i))
  preds <- unlist(sapply(fold_res, function(x) x$pl))
  true_lab <- unlist(sapply(fold_res, function(x) x$tl))
  table(data.frame(cbind(preds, true_lab )))
  sum(preds==true_lab)/length(true_lab) # 81.1%
  
  # train on all
  design <- meta_df2$sex
  exp <- expr_df4[meta_df2$sample_acc,]
  fit <- lmFit(t(exp), design)
  fit <- eBayes(fit)
  tt <-topTable(fit, number=ncol(exp))
  tt$transcript <- rownames(tt)
  top_f <- tt %>% arrange(logFC) %>% head()
  top_m <- tt %>% arrange(desc(logFC)) %>% head()
  
  list_transcripts <- c(top_f$transcript, top_m$transcript)

  fit = glmnet(expr_df4[meta_df2$sample_acc,list_transcripts], meta_df2$sex, 
               family="binomial",
               alpha=1, 
               standardize=FALSE)
  my.lambda <- 0.00208
  
  mat_coef <- coef(fit, s=my.lambda) %>% as.matrix()
  nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
  coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef))
  
  xy_genes <- read_csv(sprintf("data/rnaseq/%s/03_model_in/xy_genes_rnaseq.csv", prefix))
  
  coef_df2 <- coef_df %>% left_join(xy_genes, by=c("gene"="transcript"))
  
  coef_df3 <- coef_df2 %>%
    mutate(coef=as.numeric(as.character(coef))) %>%
    arrange(coef) %>% 
    filter(!is.na(chromosome_name)) 
  
  coef_df3 %>%
    write_csv(sprintf("data/rnaseq/%s/04_model_out/rnaseq_%s_coef.csv", prefix, prefix))
  save(fit, my.lambda, file=sprintf("data/rnaseq/%s/04_model_out/rnaseq_sl_fit.RData", prefix))
  
  # ----- what is the cv accuracy ignoring study divisions? ----- #
  # f <- meta_df2 %>% filter(sex==0)
  # m <- meta_df2 %>% filter(sex==1)
  # nfolds <- 5
  # set.seed(27)
  # f.s <- f %>% sample_n(nrow(f))
  # m.s <- m %>% sample_n(nrow(m))
  # f.s$fold <- c(rep(c(c(1:nfolds), c(nfolds:1)), floor(nrow(f.s)/(2*nfolds))))
  # # works b/c exact
  # m.s$fold <- c(rep(c(c(1:nfolds), c(nfolds:1)), floor(nrow(m.s)/(2*nfolds))),
  #               c(1:( nrow(m.s)%%(2*nfolds)  )))
  # fold_res2 <- lapply(1:nfolds, function(i) train_fold(rbind(f.s, m.s), expr_df4, i))
  # preds <- unlist(sapply(fold_res2, function(x) x$pl))
  # true_lab <- unlist(sapply(fold_res2, function(x) x$tl))
  # table(data.frame(cbind(preds, true_lab )))
  # sum(preds==true_lab)/length(true_lab) # 81.1%, 98.8% if ignore study
  
  # -- leave one study out -- #
  fold_df2 <- fold_df
  fold_df2$fold <- c(1:14)[factor(fold_df$study_acc)]
  fold_res3 <- lapply(1:14, function(i) train_fold(fold_df2, expr_df4, i))
  preds <- unlist(sapply(fold_res3, function(x) x$pl))
  true_lab <- unlist(sapply(fold_res3, function(x) x$tl))
  table(data.frame(cbind(preds, true_lab )))
  sum(preds==true_lab)/length(true_lab) # 93.9% if leave one study out
  save(fold_res3, file=sprintf("data/rnaseq/%s/04_model_out/fold_res.RData", prefix))
  
}


# --------- LET'S GO! ---------- #

# count and shuffle within each count
set.seed(10)
study_counts_df <- meta_df %>% 
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

ngroups <- 8

# need to assign 1:8 *but* also keep things fairly evenly sized
group_l <- c(1:ngroups, ngroups:1)
reps <- nrow(study_counts)%/% length(group_l)
rems <- nrow(study_counts)%% length(group_l)
opp_group_l <- c(ngroups:1, 1:ngroups)
study_counts$fold <- c(rep(group_l, reps), opp_group_l[1:rems])

# now count the number of samples, studies per group
# - this is pretty even 
study_counts %>% 
  group_by(fold) %>% 
  summarize(tot_samples=sum(num_samples), studies=n())

# sample to fold vector
samp_to_fold <- left_join(meta_df, study_counts , 
                          by=c("study_acc")) %>% select(-num_samples)

# they're fairly even, I am ok with this
samp_to_fold %>% group_by(fold) %>% 
  summarize(num_f=sum(sex=="female"), 
            num_m=sum(sex=="male"))
samp_to_fold2 <- samp_to_fold %>% 
  mutate(sex=ifelse(sex=="female", 0, 1))
test_folds <- sample(1:ngroups, 2)
test_data <- samp_to_fold2 %>% 
  filter(fold %in% test_folds)
test_expr_data <- expr_df2[,test_data$sample_acc]
test_labels <- test_data$sex

train_valid <- samp_to_fold2 %>% filter(!fold %in% test_folds)

run_fold <- function(my.fold){
  train_data <- train_valid %>% filter(!fold %in% my.fold)
  valid_data <- train_valid %>% filter(fold %in% my.fold)
  # redo the training fold
  my.folds <- data.frame(cbind("fold"=unique(train_data$fold), "new_fold"=1:5))
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
    #preds_valid <- predict(cvfit, newx=x_valid, s="lambda.1se", type="response")
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

fold_res <- data.frame()
folds <- unique(train_valid$fold)
for (i in 1:length(folds)){
  print(i)
  res <- run_fold(folds[i])
  if (i==1){
    fold_res <- res
  }
  else {
    fold_res <- rbind(fold_res, res)
  }
} 

fold_df <- apply(fold_res, c(1,2), unlist) %>% as_tibble()  %>%
  mutate(t=as.numeric(t),
         v=as.numeric(v),
         lambda=as.numeric(lambda),
         alpha=as.numeric(alpha)) 
fold_df %>% 
  write_csv(sprintf("data/rnaseq/%s/04_model_out/%s_fold_res_transform.csv", prefix, prefix))

# save the training and testing data + cutoff
train_valid_expr <- expr_df2[,train_valid$sample_acc]
train_valid_lab <- train_valid$sex
save(train_valid_expr, train_valid_lab, test_expr_data, test_labels, 
     file=sprintf("data/rnaseq/%s/03_model_in/train_test.RData", prefix))
