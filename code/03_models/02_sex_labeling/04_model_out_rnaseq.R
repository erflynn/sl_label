require('tidyverse')
require('glmnet')
require('bestNormalize')


fold_df <- read_csv(sprintf("data/rnaseq/%s/04_model_out/%s_fold_res_transform.csv", prefix, prefix))
ggplot(fold_df, aes(y=v, x=alpha, group=alpha))+
  geom_boxplot()+
  ylab("validation error")
mean_sd <- fold_df %>% 
  group_by(alpha) %>% 
  summarize(mean=mean(v), sd=sd(v)) %>%
  mutate(vlo=mean-1.96*sd, vup=mean+1.96*sd)

my.alpha <- 0.5 #0.1 for mouse 
lambda_best <- fold_df %>% 
  filter(alpha==my.alpha) %>% 
  summarize(mean=mean(lambda), sd=sd(lambda)) # --> use this!

load(sprintf("data/rnaseq/%s/03_model_in/train_test.RData", prefix))
x_train <- apply(train_valid_expr, 1, function(row) boxcox(row+0.5)$x.t)
x_test <- apply(test_expr_data, 1, function(row) boxcox(row+0.5)$x.t)

my.lambda <- lambda_best$mean
fit <- glmnet(x_train, alpha=my.alpha, train_valid_lab, family="binomial", standardize=FALSE)
train_pred <- predict(fit, newx=x_train, s=my.lambda, type="class")
sum(train_pred==train_valid_lab)/length(train_valid_lab) 

test_pred <- predict(fit, newx=x_test, s=my.lambda, type="class")
sum(test_pred==test_labels)/length(test_labels) # 0.915

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

# human only
load("../sex_labeling/geo_pipeline/gpl_ref/human_gene_map.RData")
coef_df3 %>% 
  left_join(ref_dat %>% select(ensembl_gene_id, hgnc_symbol) %>% unique()) %>% 
  filter(!is.na(hgnc_symbol)) %>% write_csv(sprintf("data/rnaseq/%s/04_model_out/coef_annot.csv", prefix))



