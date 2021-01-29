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

cvfit <- cv.glmnet(t(expr_df2[,train_filt$sample_acc]), train_filt$sex, alpha=0.5, type.measure="deviance", family="binomial")
preds_class_train <- sapply(predict(cvfit, newx=t(expr_df2[,train_filt$sample_acc]), 
                                    s="lambda.1se", type="class"), as.numeric)
sum(preds_class_train==train_filt$sex)/length(train_filt$sex) 
preds_train <- predict(cvfit, newx=t(expr_df2[,train_filt$sample_acc]), s="lambda.1se", type="response")

preds_test <- predict(cvfit, newx=t(expr_df2[,test_filt$sample_acc]), s="lambda.1se", type="response")
preds_class_test <- sapply(predict(cvfit, newx=t(expr_df2[,test_filt$sample_acc]), 
                                   s="lambda.1se", type="class"), as.numeric)
sum(preds_class_test==test_filt$sex)/length(test_filt$sex) 
plot(density(preds_train[preds_class_train != train_filt$sex,]))

plot(density(preds_test[preds_class_test != test_filt$sex,]))


read_counts_sm <- read_counts %>% filter(sample_acc %in% c(colnames(train_valid_expr), colnames(test_expr_data)))
