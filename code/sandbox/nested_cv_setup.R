# setup:
# -- meta_df: study, sample, label (binary for now)
# -- expr_df: cols are samples
# goal:
# each study is in a unique fold
#  80% for nested CV
#   8-fold CV
#    1-6 --> into cv.glmnet to tune for lambda, alpha, type.measure, etc
#    valid on two folds (before too small)
#  20% held out test

require('tidyverse')
require('glmnet')

## ------- SETUP -------- ##

load(sprintf("data/rnaseq/%s/03_model_in/expr_sl_in.RData", prefix)) # --> all_df

sl <- read_csv(sprintf("data/01_metadata/%s_rnaseq_sex_lab.csv",prefix)) %>% 
  as.data.frame()
rownames(sl) <- sl$sample_acc

# filter so no cell line samples (?)
metadata <- read.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", prefix), stringsAsFactors = FALSE)
metadata_sm <- metadata %>% filter(acc %in% colnames(all_df))
no_cl <- metadata_sm %>% 
  filter(cl_line %in% c("", "--", "n/a", "not applicable", "na"))

samples <- intersect(no_cl$acc, sl$sample_acc)
expr_df <- all_df[,samples]
meta_df <- sl[samples,] 

# setup genes
xy_genes <- read_csv(sprintf("data/rnaseq/%s/03_model_in/xy_genes_rnaseq.csv", prefix))
rownames(expr_df) <- all_df$gene_name
expr_df2 <- expr_df[xy_genes$transcript,]

# --------- LET'S GO! ---------- #

# count and shuffle within each count
set.seed(5)
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

# need to assign 1:8 *but* also keep things fairly evenly sized
group_l <- c(1:8, 8:1)
reps <- nrow(study_counts)%/% length(group_l)
rems <- nrow(study_counts)%% length(group_l)
opp_group_l <- c(8:1, 1:8)
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
test_folds <- sample(1:8, 2)
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
  
  # now run the fit!
  x_train <- t(train_expr_data)
  x_valid <- t(valid_expr_data)

  test_params <- function(my.alpha, my.measure){
    cvfit = cv.glmnet(x_train, train_labels, 
                      family="binomial",
                      foldid=train_folds,
                      nfolds=5,
                      alpha=my.alpha, 
                      standardize=TRUE,
                      type.measure=my.measure)
    preds_train <- predict(cvfit, newx=x_train, s="lambda.1se", type="response")
    preds_class_train <- sapply(predict(cvfit, newx=x_train, s="lambda.1se", type="class"), 
                                as.numeric)
    train_acc <- sum(preds_class_train==train_labels)/length(train_labels)
    #preds_valid <- predict(cvfit, newx=x_valid, s="lambda.1se", type="response")
    preds_class_valid <- sapply(predict(cvfit, newx=x_valid, s="lambda.1se", type="class"), as.numeric)
    valid_acc <- sum(preds_class_valid==valid_labels)/length(valid_labels) 
    return(list("t"=train_acc, 
                "v"=valid_acc, 
                "lambda"=cvfit$lambda.1se, 
                "alpha"=my.alpha, 
                "measure"=my.measure))
  }
  param_res <- data.frame()
  for (my.alpha in seq(0,1, 0.1)){
    alpha_res <- test_params(my.alpha, "deviance")
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
  write_csv(sprintf("data/rnaseq/%s/04_model_out/%s_fold_res.csv", prefix, prefix))

# save the training and testing data + cutoff
train_valid_expr <- expr_df2[,train_valid$sample_acc]
train_valid_lab <- train_valid$sex
save(train_valid_expr, train_valid_lab, test_expr_data, test_labels, 
     file=sprintf("data/rnaseq/%s/03_model_in/train_test.RData", prefix))


# setup the model in a function

ggplot(fold_df, aes(y=v, x=alpha, group=alpha))+
  geom_boxplot(aes(col=measure))+
  ylab("validation error")
mean_sd <- fold_df %>% 
  group_by(alpha, measure) %>% 
  summarize(mean=mean(v), sd=sd(v)) %>%
  mutate(vlo=mean-1.96*sd, vup=mean+1.96*sd)
mean_sd %>% filter(measure=="deviance")
fold_df %>% 
  filter(alpha==0.1, measure=="deviance") %>% 
  summarize(mean=mean(lambda), sd=sd(lambda)) # --> use this!

# alpha=0.3, measure="deviance"
my.folds <- data.frame(cbind("fold"=unique(train_valid$fold), "new_fold"=1:6))
train_data2 <- train_valid %>% left_join(my.folds)
train_folds <- train_data2$new_fold

cvfit <- cv.glmnet(t(train_valid_expr), alpha=1, foldid=train_folds, train_valid_lab, family="binomial")
train_pred <- predict(cvfit, newx=t(train_valid_expr), s="lambda.1se", type="class")
sum(train_pred==train_valid_lab)/length(train_valid_lab) 

test_pred <- predict(cvfit, newx=t(test_expr_data), s="lambda.1se", type="class")
sum(test_pred==test_labels)/length(test_labels) # 0.915

mat_coef <- coef(cvfit, lambda="lambda.1se") %>% as.matrix()
nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef))

xy_genes <- read_csv(sprintf("data/rnaseq/%s/03_model_in/xy_genes_rnaseq.csv", prefix))

coef_df2 <- coef_df %>% left_join(xy_genes, by=c("gene"="transcript"))

coef_df2 %>%
  mutate(coef=as.numeric(as.character(coef))) %>%
  arrange(coef) %>% 
  filter(!is.na(chromosome_name)) %>% 
  write_csv(sprintf("data/%s/04_sl_input/%s_coef.csv", prefix, prefix))

save(cvfit, file=sprintf("data/%s/04_sl_input/cvfit.RData", prefix))

coef_df2 %>% arrange(coef) %>% left_join(ref_dat %>% 
                                          select(ensembl_gene_id, hgnc_symbol) %>% 
                                          unique(), by=c("gene"="ensembl_gene_id")) %>% 
  filter(!is.na(hgnc_symbol)) %>% write_csv(sprintf("data/%s/04_sl_input/coef_annot.csv", prefix))


save(fit, file=sprintf("data/rnaseq/%s/04_model_out/%s_rnaseq_fit.RData", prefix, prefix))



