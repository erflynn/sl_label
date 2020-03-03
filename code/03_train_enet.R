require('glmnet')
require('tidyverse')
require('data.table')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

sex_lab_f <- fread(sprintf("data/%s/04_sl_input/%s_f_sex_lab.csv", prefix, prefix), 
                   data.table=FALSE, header=TRUE)
sex_lab_m <- fread(sprintf("data/%s/04_sl_input/%s_m_sex_lab.csv", prefix, prefix ), 
                   data.table=FALSE, header=TRUE)
sex_lab_f$V1 <- NULL
sex_lab_m$V1 <- NULL

expr_f <- fread(sprintf("data/%s/04_sl_input/%s_f_expr.csv", prefix, prefix ), data.table=FALSE)
expr_m <- fread(sprintf("data/%s/04_sl_input/%s_m_expr.csv", prefix, prefix ), data.table=FALSE)

table(colnames(expr_f[,2:1001])==sex_lab_f$acc)
table(colnames(expr_m[,2:1001])==sex_lab_m$acc)

# make sure the data are ok
load(sprintf("../sex_labeling/geo_pipeline/gpl_ref/%s_gene_map.RData", prefix))
xy_genes <- ref_dat %>% 
  filter(chromosome_name %in% c("X", "Y")) %>% 
  select(chromosome_name, ensembl_gene_id) %>% 
  unique()

gene_names <- expr_f$rid
overlap.genes <- intersect(gene_names,xy_genes$ensembl_gene_id)
rownames(expr_f) <- gene_names
expr_f <- expr_f[overlap.genes,]
rownames(expr_m) <- gene_names
expr_m <- expr_m[overlap.genes,]

# train/test
set.seed(301)
sample_f <- sample(1:1000, 700, replace=FALSE)
sample_m <- sample(1:1000, 700, replace=FALSE)

expr_f_train <- (expr_f %>% select(-rid))[,sample_f]
expr_m_train <- (expr_m %>% select(-rid))[,sample_m]

expr_f_test <- (expr_f %>% select(-rid))[,-sample_f]
expr_m_test <- (expr_m %>% select(-rid))[,-sample_m]


expr_train <- cbind(expr_f_train, expr_m_train)
train_sex_lab <- c(rep(0, 700), rep(1, 700))

expr_test <- cbind(expr_f_test, expr_m_test)
test_sex_lab <- c(rep(0, 300), rep(1, 300))

# train lasso on these data
x_train <- as.matrix(t(expr_train))
x_train <- apply(x_train, c(1,2),as.numeric)

x_test <- as.matrix(t(expr_test))
x_test <- apply(x_test, c(1,2),as.numeric)

cvfit = cv.glmnet(x_train, train_sex_lab, family="binomial", alpha=0.5)
preds_train <- predict(cvfit, newx=x_train, s="lambda.min", type="response")
preds_class_train <- sapply(predict(cvfit, newx=x_train, s="lambda.min", type="class"), as.numeric)
sum(preds_class_train==train_sex_lab)/length(train_sex_lab) # 96.4%

preds_test <- predict(cvfit, newx=x_test, s="lambda.min", type="response")
preds_class_test <- sapply(predict(cvfit, newx=x_test, s="lambda.min", type="class"), as.numeric)
sum(preds_class_test==test_sex_lab)/length(test_sex_lab) # 95.5%

mat_coef <- coef(cvfit, lambda="lambda.min") %>% as.matrix()
nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef))
coef_df2 <- coef_df %>% left_join(xy_genes, by=c("gene"="ensembl_gene_id"))

coef_df2 %>% 
  arrange(coef) %>% 
  filter(!is.na(chromosome_name)) %>% 
  write_csv(sprintf("data/%s/04_sl_input/%s_coef.csv", prefix, prefix))

save(cvfit, file=sprintf("data/%s/04_sl_input/cvfit.RData", prefix))

coef_df2 %>% arrange(coef) %>% left_join(ref_dat %>% 
                                          select(ensembl_gene_id, hgnc_symbol) %>% 
                                          unique(), by=c("gene"="ensembl_gene_id")) %>% 
  filter(!is.na(hgnc_symbol)) %>% write_csv(sprintf("data/%s/04_sl_input/coef_annot.csv", prefix))
