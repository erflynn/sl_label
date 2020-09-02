
require('bestNormalize')
miceadds::load.Rdata(sprintf("data/05_train_df/%s_rnaseq_%s_pos_expr.RData", prefix, ds), "pos")
miceadds::load.Rdata(sprintf("data/05_train_df/%s_rnaseq_%s_neg_expr.RData", prefix, ds), "neg")
pos <- pos %>% rename(rid=gene_name)
neg <- neg %>% rename(rid=gene_name) 
load(file=sprintf("data/07_model_dat/fit_%s_rnaseq_%s.RData", prefix, ds))
list_genes <- unique(rownames(fit$beta))
dat <- full_join(pos, neg, by="rid") %>% as_tibble() %>% filter(rid %in% list_genes)


load(sprintf("data/07_model_dat/%s_rnaseq_sex_train_mat.RData", prefix))
train_data <- rownames(X_train)
test_data <- rownames(X_test)

train <- dat[,train_data]
test <- dat[,test_data]
rm(X_train, X_test)
X_train <- t(train)
X_test <- t(test)
colnames(X_train) <- dat$rid
colnames(X_test) <- dat$rid
X_train[is.na(X_train)] <- 0
X_test[is.na(X_test)] <- 0
box_cox_obj <- apply(X_train, 2, function(col) boxcox(col+0.5)) 
X_test2.1 <- predict(box_cox_obj[[1]], newx=(X_test[,1]+0.5))
X_train2 <- do.call(cbind, lapply(box_cox_obj, function(x) x$x.t) )
X_test2 <- do.call(cbind, lapply(1:length(list_genes), 
                                 function(i) predict(box_cox_obj[[i]], 
                                                     newdata=(X_test[,i]+0.5) )
))
#X_test <- apply(X_test, 2, function(col) boxcox(col+0.5)$x.t)

save(list_genes, box_cox_obj, 
     file=sprintf("data/07_model_dat/%s_rnaseq_boxcox.RData", prefix))