library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
data_type <- args[2]
#ds <- args[3]

load(sprintf("data/data_old/07_model_dat/%s_%s_sex_train_mat.RData", prefix, data_type))

sex_lab <- ifelse(Y_train==0, "female", "male")

# X_test, X_train, Y_test, Y_train
#train_df <- cbind(X_train, "sex"=ifelse(Y_train==0, "female", "male"))
# https://machinelearningmastery.com/machine-learning-in-r-step-by-step/

library(caret)
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"

# LDA
print("LDA")
set.seed(1)
fit.lda <- train(X_train, sex_lab, method="lda", metric=metric, 
                 trControl=control)
# CART
print("CART")
set.seed(1)
fit.cart <- train(X_train, sex_lab, method="rpart", metric=metric, 
                  trControl=control)
# kNN
print("KNN")
set.seed(1)
fit.knn <- train(X_train, sex_lab, method="knn", metric=metric, 
                 trControl=control)

print("SVM")
set.seed(1)
fit.svm <- train(X_train, sex_lab, method="svmRadial", metric=metric, 
                 trControl=control)
# Random Forest
print("RF")
set.seed(1)
fit.rf <- train(X_train, sex_lab, 
                method="rf", metric=metric, 
                trControl=control)

library('randomForest')
rf_classifier = randomForest(x=X_train, y=factor(sex_lab), ntree=500)
pred = predict(rf_classifier, X_test)
sum(ifelse(pred=="female", 0, 1)==Y_test)/length(Y_test) # 92.5%
  
# GLM
print("Glmnet")
set.seed(1)
#glmGrid <- expand.grid(alpha=c(0,0.5, 1), lambda=0.05)
fit.glmnet <-  train(X_train, sex_lab, 
                  method="glmnet", metric=metric,
                  trControl=control)#,
#                  tuneGrid=glmGrid)


# method="glm"
# TODO: update for mouse, RNA-seq 
####
print("LM")
library(limma)
sex = factor(Y_train)
design_l <- model.matrix(~sex)
fit_l <- lmFit(t(X_train), design_l)
fit_l <- eBayes(fit_l)
tt <- data.frame(topTable(fit_l, n=20))
tt$ensembl_gene_id <- rownames(tt) 
topy <- tt %>% arrange(desc(logFC)) %>% head(1) %>% pull(ensembl_gene_id)
topx <- tt %>% arrange(desc(-logFC)) %>% head(1) %>% pull(ensembl_gene_id)
topx
topy
train_df2 <- data.frame(train_df[,c(topy, topx, "sex")]) %>%
  mutate(across(c(topy, topx), ~as.numeric(as.character(.))))
fit.lm <- train(train_df2[,c(1:2)], train_df2$sex, 
                method="glm", metric=metric,
                trControl=control)

save(fit.cart, fit.lda, fit.knn, fit.svm, fit.glmnet, fit.rf, fit.lm,
     file=sprintf("data/10_ml_models/%s_%s_alt_ml_models.RData", prefix, data_type))
results <- resamples(list(lda=fit.lda, cart=fit.cart, 
                          knn=fit.knn, svm=fit.svm, rf=fit.rf,
                          lm=fit.lm, glmnet=fit.glmnet))
print(summary(results))


load("data/alt_ml_models.RData")
library(caret)
results <- resamples(list(lda=fit.lda, cart=fit.cart, 
                          knn=fit.knn, svm=fit.svm,
                          lm=fit.lm, glmnet=fit.glmnet))
print(summary(results))
#### other setup

cvfit <- cv.glmnet(X_train, Y_train, family="binomial", alpha=0.5, nfolds=6)
summary(cvfit)
train_assess <- assess.glmnet(cvfit, newx=X_train, newy=Y_train)
test_pred <- predict(cvfit, newx=X_test, type="response")
assess.glmnet(cvfit, newx=X_test, newy=Y_test)# 
# # plot AUC?
# 
# # --- TODO: wrap everything in CV --- # 
# 
# # PC regression
# pcs <- prcomp(X_train)
# plot(pcs$x[,"PC1"], pcs$x[,"PC2"], 
#      col=c("red", "blue")[Y_train])
# 
# # PC regression w/ X + Y
# # PC regression w/ X + Y divided
# 
# 
# 
# # logistic regression with top coef
# 
# load("../cersi_code/data/rb_ensembl_to_hgnc.RData")
# tt %>% left_join(convert_genes, by="ensembl_gene_id")
# 
# # RPS4Y1, KDM5D, XIST
# fit_s <- glm(factor(sex) ~ ENSG00000129824 + ENSG00000012817+ ENSG00000229807 ,
#              data=data.frame(train_df), family="binomial")
# summary(fit_s)
# preds <- predict(fit_s, newx=X_train, type="response")
# preds[preds<0.5] <- 0
# preds_c <- ifelse(preds < 0.5, 0, 1)
# table(preds_c==Y_train)
# 
# sum(preds_c==Y_train)/length(preds_c) # 89.8
# 
# 
# 
# fit_s <- glm(factor(sex) ~ ENSG00000129824 + ENSG00000229807,
#              data=data.frame(train_df), family="binomial")
# summary(fit_s)
# preds <- predict(fit_s, newx=X_train, type="response")
# preds[preds<0.5] <- 0
# preds_c <- ifelse(preds < 0.5, 0, 1)
# table(preds_c==Y_train)
# 
# sum(preds_c==Y_train)/length(preds_c) # 84.6
# 
# # glmnet w/ alpha = 0, 0.5, 1 w/o cv
# library(glmnet)
# cvfit = cv.glmnet(X_train, Y_train, 
#                   family="binomial",
#                   #foldid=train_folds,
#                   #alpha=my.alpha, 
#                   alpha=0.5,
#                   standardize=FALSE,
#                   #standardize=standardizeFlag,
#                   trace=FALSE)
# (train_assess <- assess.glmnet(cvfit, newx=X_train, newy=Y_train))
# 
# cvfit$lambda
# 
# 
# cvfit0 = cv.glmnet(X_train, Y_train, 
#                   family="binomial",
#                   #foldid=train_folds,
#                   #alpha=my.alpha, 
#                   alpha=0,
#                   standardize=FALSE,
#                   #standardize=standardizeFlag,
#                   trace=FALSE)
# (train_assess <- assess.glmnet(cvfit0, newx=X_train, newy=Y_train))
# 
# cvfit1 = cv.glmnet(X_train, Y_train, 
#                   family="binomial",
#                   #foldid=train_folds,
#                   #alpha=my.alpha, 
#                   alpha=1,
#                   standardize=FALSE,
#                   #standardize=standardizeFlag,
#                   trace=FALSE)
# (train_assess <- assess.glmnet(cvfit1, newx=X_train, newy=Y_train))
# 
# 
# valid_pred <- predict(cvfit, newx=X_valid2, type="response")
# valid_assess <- assess.glmnet(cvfit, newx=X_valid2, newy=Y_valid2)
# 
# 
# 
# # estimate skill of LDA on the validation dataset
# #validation_df <- data.frame(cbind(X_test, "sex"=Y_test))
# 
# #predictions <- predict(fit.cart, X_test)
# #confusionMatrix(predictions, factor(validation_df$sex))
# 
# 
# 
# # other methods:
# #  - deep learning
