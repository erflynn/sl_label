require('glmnet')
require('tidyverse')
require('data.table')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

expr_cl <- fread(sprintf("data/%s/04_sl_input/cell_pos_expr.csv", prefix, prefix ), data.table=FALSE)
expr_no_cl <- fread(sprintf("data/%s/04_sl_input/cell_neg_expr.csv", prefix, prefix ), data.table=FALSE)

expr_cl_train <- expr_cl[,2:3000]
expr_cl_test <- expr_cl[,3001:6000]

expr_no_cl_train <- expr_no_cl[,2:3000]
expr_no_cl_test <- expr_no_cl[,3001:6000]

expr_train <- cbind(expr_no_cl_train, expr_cl_train)
rownames(expr_train) <- expr_no_cl$rid

train_cl_lab <- c(rep(0, ncol(expr_no_cl_train)), rep(1, ncol(expr_cl_train)))

expr_test <- cbind(expr_no_cl_test, expr_cl_test)
rownames(expr_test) <- expr_no_cl$rid
test_cl_lab <- c(rep(0, ncol(expr_no_cl_test)), rep(1, ncol(expr_cl_test)))


# train lasso on these data
x_train <- as.matrix(t(expr_train))
#x_train <- apply(x_train, c(1,2),as.numeric)

x_test <- as.matrix(t(expr_test))
#x_test <- apply(x_test, c(1,2),as.numeric)
y <-train_cl_lab
foldid=sample(1:10,size=length(y),replace=TRUE)
x <- x_train
cv1=cv.glmnet(x,y,foldid=foldid,alpha=1, family="binomial", standardize=FALSE, type.measure="class")
cv.5=cv.glmnet(x,y,foldid=foldid,alpha=.5, family="binomial", standardize=FALSE, type.measure="class")
cv0=cv.glmnet(x,y,foldid=foldid,alpha=0, family="binomial", standardize=FALSE, type.measure="class")
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
