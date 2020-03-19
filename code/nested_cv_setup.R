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

exp_sample <- read.csv("data/01_metadata/human_rnaseq_exp_to_sample2.csv")
metadata <- read.csv("data/01_metadata/human_rnaseq_sample_metadata.csv")
load("human_xy_expr_df.RData")
prefix <- "human"
sl <- read_csv(sprintf("data/01_metadata/%s_rnaseq_sex_lab.csv",  prefix))
sl_sm <- sl %>% filter(sample_acc %in% colnames(sample_expr_df2)) %>% as.data.frame()
metadata_sm <- metadata %>% filter(acc %in% colnames(sample_expr_df2))
no_cl <- metadata_sm %>% filter(cl_line=="")
keep.cols <- colnames(sample_expr_df2) %in% no_cl$acc
sample_expr_df3 <- sample_expr_df2[,keep.cols]
rownames(sl_sm) <- sl_sm$sample_acc
sl_sm2 <- sl_sm[colnames(sample_expr_df3),]

meta_df <- sl_sm2
expr_df <- sample_expr_df3

# --------- LET'S GO! ---------- #

# count and shuffle within each count
set.seed(5)
study_counts_df <- meta_df %>% group_by(study_acc) %>% count() %>% arrange(desc(n)) %>%
  ungroup() %>% rename(num_samples=n) %>% group_by(num_samples) %>% mutate(n=n() ) %>%
  sample_n(n) %>% mutate(ord=1:n()) %>% arrange(desc(num_samples, ord))

study_counts <- study_counts_df %>% select(study_acc, num_samples) %>% ungroup()

# need to assign 1:10 *but* also keep things fairly evenly sized
group_l <- c(1:10, 10:1)
reps <- nrow(study_counts_df)%/% length(group_l)
rems <- nrow(study_counts_df)%% length(group_l)
opp_group_l <- c(10:1, 1:10)
study_counts$fold <- c(rep(group_l, reps), opp_group_l[1:rems])

# now count the number of samples, studies per group
# - this is pretty even 
study_counts %>% group_by(fold) %>% summarize(tot_samples=sum(num_samples), studies=n())

# sample to fold vector
samp_to_fold <- left_join(meta_df, study_counts, by=c("study_acc")) %>% select(-num_samples)

# they're fairly even, I am ok with this
samp_to_fold %>% group_by(fold) %>% summarize(num_f=sum(sex=="female"), num_m=sum(sex=="male"))
samp_to_fold2 <- samp_to_fold %>% mutate(sex=ifelse(sex=="female", 0, 1))
test_folds <- sample(1:10, 2)
test_data <- samp_to_fold2 %>% filter(fold %in% test_folds)
test_expr_data <- expr_df[,test_data$sample_acc]
test_labels <- test_data$sex
  
train_valid <- samp_to_fold2 %>% filter(!fold %in% test_folds)

run_fold <- function(my.fold){
  train_data <- train_valid %>% filter(!fold %in% my.fold)
  valid_data <- train_valid %>% filter(fold %in% my.fold)
  # redo the training fold
  my.folds <- data.frame(cbind("fold"=unique(train_data$fold), "new_fold"=1:6))
  train_data2 <- train_data %>% left_join(my.folds)
  train_expr_data <- expr_df[,train_data2$sample_acc]
  train_labels <- train_data2$sex
  train_folds <- train_data2$new_fold
  valid_expr_data <- expr_df[,valid_data$sample_acc]
  valid_labels <- valid_data$sex
  
  # now run the fit!
  x_train <- t(train_expr_data)
  x_valid <- t(valid_expr_data)
  #my.alpha <- 0.5
  #my.measure <- "deviance" #c("class", "auc", "mse", "mae")
  test_params <- function(my.alpha, my.measure){
    cvfit = cv.glmnet(x_train, train_labels, 
                      family="binomial",
                      foldid=train_folds,
                      nfolds=7,
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
    r1 <- test_params(my.alpha, "class")
    r2 <- test_params(my.alpha, "deviance")
    alpha_res <- rbind(r1, r2)
    if (my.alpha==0){
      param_res <- alpha_res
    } else {
      param_res <- rbind(param_res, alpha_res)
    }
  }
  return(param_res)
}

combs <- combn(unique(train_valid$fold),2)
for (i in 1:ncol(combs)){
  print(i)
  res <- run_fold(combs[,i])
  if (i==1){
    fold_res <- res
  }
  else {
    fold_res <- rbind(fold_res, res)
  }
} 

fold_res <- lapply(unique(train_valid$fold), run_fold)
fold_df <- apply(fold_res, c(1,2), unlist) %>% as_tibble()  %>%
  mutate(t=as.numeric(t),
         v=as.numeric(v),
         lambda=as.numeric(lambda),
         alpha=as.numeric(alpha),
         measure=as.factor(measure)) %>%
  mutate(my_grp=paste(measure,alpha, sep="."))
ggplot(fold_df, aes(y=v, x=alpha, group=my_grp))+
  #geom_violin(aes(col=measure), draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_boxplot(aes(col=measure))+
  ylab("validation error")
mean_sd <- fold_df %>% 
  group_by(alpha, measure) %>% 
  summarize(mean=mean(v), sd=sd(v)) %>%
  mutate(vlo=mean-1.96*sd, vup=mean+1.96*sd)
mean_sd %>% filter(measure=="deviance")
fold_df %>% filter(alpha==0.3, measure=="deviance") %>% summarize(mean=mean(lambda), sd=sd(lambda))

# alpha=0.3, measure="deviance"
train_valid_expr <- expr_df[,train_valid$sample_acc]
train_valid_lab <- train_valid$sex
fit <- glmnet(t(train_valid_expr), alpha=0.3, train_valid_lab, family="binomial")
mat_coef <- coef(fit,s=0.130)  %>% as.matrix()
nonzero_coef <- mat_coef[mat_coef[,1]!=0,] # 61
train_valid_pred <- predict(fit, newx=t(train_valid_expr), s=0.130, type="class")
sum(train_valid_pred==train_valid_lab)/length(train_valid_lab) # 0.863

test_pred <- predict(fit, newx=t(test_expr_data), s=0.130, type="class")
sum(test_pred==test_labels)/length(test_labels) # 0.939

test_pv <- predict(fit, newx=t(test_expr_data), s=0.130, type="response") 
test_pred_val <- data.frame("pred"=test_pv, "acc"=rownames(test_pv))
test_pred_val$true_lab <- test_labels
train_pv <- predict(fit, newx=t(train_valid_expr), s=0.130, type="response")
train_pred_val <- data.frame("pred"=train_pv, "acc"=rownames(train_pv))

train_pred_val$true_lab <- train_valid_lab
pred_df <- rbind(train_pred_val %>% mutate(src="train"), test_pred_val %>% 
                   mutate(src="test"))%>% rename(pred=X1) %>% mutate(pred_lab=ifelse(pred < 0.5, 0, 1)) %>%
  as_tibble()




# LOL...
# --- bigger held out test?
read_counts <- read_csv("data/all_sl_read_counts.csv")

pred_df2 <- pred_df %>% left_join(read_counts, by=c("acc"="sample_acc")) %>% mutate(correct=(true_lab==pred_lab))
pred_df2 %>% group_by(correct) %>% summarize(my_min=min(num_reads), 
                                             median=median(num_reads),
                                             per1=quantile(num_reads, probs=0.01),
                                             per5=quantile(num_reads, probs=0.05),
                                             per10=quantile(num_reads, probs=0.1))
pred_df2 %>% summarize(my_min=min(num_reads), 
                                             median=median(num_reads),
                                             per1=quantile(num_reads, probs=0.01),
                                             per5=quantile(num_reads, probs=0.05),
                                             per10=quantile(num_reads, probs=0.1))
ggplot(pred_df2 %>% filter(src=="train"), aes(x=correct, y=num_reads))+
  geom_violin(aes(col=correct), draw_quantiles=c(0.25, 0.5, 1)) 

(pred_df2 %>% filter(num_reads >= 100000) %>% nrow())/(pred_df2 %>% nrow()) # 97.2%
ggplot(pred_df2 %>% filter(num_reads >= 100000), aes(x=abs(0.5-pred), y=num_reads))+
  geom_point(aes(col=correct), alpha=0.5)

ggplot(pred_df2 %>% filter(num_reads >= 100000), aes(x=abs(0.5-pred)))+
  geom_density(aes(col=correct))
pred_df2 %>% filter(num_reads >=500000) %>% filter(src=="train") %>% summarize(sum(correct)/n())
pred_df2 %>% filter(num_reads >=500000) %>% filter(src=="test") %>% summarize(sum(correct)/n())
samp_to_fold3 <- samp_to_fold2 %>% left_join(read_counts %>% select(sample_acc, num_reads)) %>%
  filter(num_reads >=100000)
set.seed(6)
test.folds <- sample(1:10, 5)
train_filt <- samp_to_fold3 %>% filter(!fold %in% test.folds)
test_filt <- samp_to_fold3 %>% filter(fold %in% test.folds)

cvfit <- cv.glmnet(t(expr_df[,train_filt$sample_acc]), train_filt$sex, alpha=0.5, type.measure="deviance", family="binomial")
preds_class_train <- sapply(predict(cvfit, newx=t(expr_df[,train_filt$sample_acc]), 
                                    s="lambda.1se", type="class"), as.numeric)
sum(preds_class_train==train_filt$sex)/length(train_filt$sex) 
preds_train <- predict(cvfit, newx=t(expr_df[,train_filt$sample_acc]), s="lambda.1se", type="response")

preds_test <- predict(cvfit, newx=t(expr_df[,test_filt$sample_acc]), s="lambda.1se", type="response")
preds_class_test <- sapply(predict(cvfit, newx=t(expr_df[,test_filt$sample_acc]), 
                                   s="lambda.1se", type="class"), as.numeric)
sum(preds_class_test==test_filt$sex)/length(test_filt$sex) 
plot(density(preds_train[preds_class_train != train_filt$sex,]))

plot(density(preds_test[preds_class_test != test_filt$sex,]))


read_counts_sm <- read_counts %>% filter(sample_acc %in% c(colnames(train_valid_expr), colnames(test_expr_data)))
