
# tried this but not sure it's right:
#  https://rpubs.com/kaz_yos/alasso
require('glmnet')
data(BinomialExample)

# adaptive
x_bin <- x
y_bin <- y
## Perform ridge regression
ridge3 <- glmnet(x = x_bin, y = y_bin,
                 ## Binary logistic regression
                 family = "binomial",
                 ## ‘alpha = 1’ is the lasso penalty, and ‘alpha = 0’ the ridge penalty.
                 alpha = 0)
plot(ridge3, xvar = "lambda")
ridge3_cv <- cv.glmnet(x = x_bin, y = y_bin,
                       ## type.measure: loss to use for cross-validation.
                       type.measure = "deviance",
                       ## K = 10 is the default.
                       nfold = 10,
                       ## Multinomial regression
                       family = "binomial",
                       ## ‘alpha = 1’ is the lasso penalty, and ‘alpha = 0’ the ridge penalty.
                       alpha = 0)
## Penalty vs CV MSE plot
plot(ridge3_cv)
(best_ridge_coef3 <- coef(ridge3_cv, s = ridge3_cv$lambda.min))
best_ridge_coef3 <- as.numeric(best_ridge_coef3)[-1] # remove intercept
## Perform adaptive LASSO
alasso3 <- glmnet(x = x_bin, y = y_bin,
                  ## Multinomial regression
                  family = "binomial",
                  ## ‘alpha = 1’ is the lasso penalty, and ‘alpha = 0’ the ridge penalty.
                  alpha = 1,
                  ##
                  ## penalty.factor: Separate penalty factors can be applied to each
                  ##           coefficient. This is a number that multiplies ‘lambda’ to
                  ##           allow differential shrinkage. Can be 0 for some variables,
                  ##           which implies no shrinkage, and that variable is always
                  ##           included in the model. Default is 1 for all variables (and
                  ##           implicitly infinity for variables listed in ‘exclude’). Note:
                  ##           the penalty factors are internally rescaled to sum to nvars,
                  ##           and the lambda sequence will reflect this change.
                  penalty.factor = 1 / abs(best_ridge_coef3))
plot(alasso3, xvar = "lambda")

lasso3 <- glmnet(x = x_bin, y = y_bin,
                 ## Multinomial regression
                 family = "binomial",
                 ## ‘alpha = 1’ is the lasso penalty, and ‘alpha = 0’ the ridge penalty.
                 alpha = 1)

lasso3_cv <- cv.glmnet(x = x_bin, y = y_bin,
                       ## type.measure: loss to use for cross-validation.
                       type.measure = "deviance",
                       ## K = 10 is the default.
                       nfold = 10,
                       ## Multinomial regression
                       family = "binomial",
                       ## ‘alpha = 1’ is the lasso penalty, and ‘alpha = 0’ the ridge penalty.
                       alpha = 1, 
                       keep=TRUE)


## Perform adaptive LASSO with 10-fold CV
alasso3_cv <- cv.glmnet(x = x_bin, y = y_bin,
                        ## type.measure: loss to use for cross-validation.
                        type.measure = "deviance",
                        ## K = 10 is the default.
                        nfold = 10,
                        ## Multinomial regression
                        family = "binomial",
                        ## ‘alpha = 1’ is the lasso penalty, and ‘alpha = 0’ the ridge penalty.
                        alpha = 1,
                        ##
                        ## penalty.factor: Separate penalty factors can be applied to each
                        ##           coefficient. This is a number that multiplies ‘lambda’ to
                        ##           allow differential shrinkage. Can be 0 for some variables,
                        ##           which implies no shrinkage, and that variable is always
                        ##           included in the model. Default is 1 for all variables (and
                        ##           implicitly infinity for variables listed in ‘exclude’). Note:
                        ##           the penalty factors are internally rescaled to sum to nvars,
                        ##           and the lambda sequence will reflect this change.
                        penalty.factor = 1 / abs(best_ridge_coef3),
                        keep = TRUE)
## Penalty vs CV MSE plot
plot(alasso3_cv)
coef(alasso3_cv, s = alasso3_cv$lambda.min)

## Extract predicted probabilities and observed outcomes.
pY <- as.numeric(predict(alasso3, newx = x_bin, s = alasso3_cv$lambda.min, type = "response"))
pY2 <- as.numeric(predict(lasso3, newx=x_bin, s=lasso3_cv$lambda.min, type="response"))
Y <- as.numeric(y_bin)
## pROC for ROC construction
roc1 <- pROC::roc(Y ~ pY)
## Plot an ROC curve with AUC and threshold
plot(roc1, print.auc = TRUE, print.thres = TRUE, print.thres.best.method = "youden")
roc2 <- pROC::roc(Y ~ pY2)
plot(roc2, print.auc = TRUE, print.thres = TRUE, print.thres.best.method = "youden")


alasso3_res <- lapply(unique(alasso3_cv$foldid), function(id) {
  ## Fit excluding test set (foldid == id)
  fit <- glmnet(x = x_bin[alasso3_cv$foldid != id,],
                y = y_bin[alasso3_cv$foldid != id],
                family = "binomial",
                alpha = 1,
                penalty.factor = 1 / abs(best_ridge_coef3))
  ## Test-set Y_hat using model fit at best lambda
  y_pred <- as.vector(predict(fit, newx = x_bin[alasso3_cv$foldid == id,], s = alasso3_cv$lambda.1se))
  ## Test-set Y
  y <- y_bin[alasso3_cv$foldid == id]
  ## Test-set AUC
  return(list("y"=y, "y_pred"=y_pred))
})
lasso3_res <- lapply(unique(lasso3_cv$foldid), function(id) {
  ## Fit excluding test set (foldid == id)
  fit <- glmnet(x = x_bin[lasso3_cv$foldid != id,],
                y = y_bin[lasso3_cv$foldid != id],
                family = "binomial",
                alpha = 1)
  ## Test-set Y_hat using model fit at best lambda
  y_pred <- as.vector(predict(fit, newx = x_bin[lasso3_cv$foldid == id,], s = lasso3_cv$lambda.1se))
  ## Test-set Y
  y <- y_bin[lasso3_cv$foldid == id]
  ## Test-set AUC
  return(list("y"=y, "y_pred"=y_pred))
}) 
all_y <- unlist(sapply(lasso3_res, function(x) x$y))
all_y_pred <- unlist(sapply(lasso3_res, function(x) x$y_pred))
pROC::roc(all_y ~ all_y_pred)$auc

all_ay <- unlist(sapply(alasso3_res, function(x) x$y))
all_ay_pred <- unlist(sapply(alasso3_res, function(x) x$y_pred))
pROC::roc(all_ay ~ all_ay_pred)$auc
