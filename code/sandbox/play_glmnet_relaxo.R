

require('glmnet')
# try BinomialExample
data(BinomialExample)
res <- cv.glmnet(x, y, family="binomial", trace=TRUE, 
                 relax=TRUE, path=TRUE)
plot(res)
# this still looks like a mess except gamma=1... hmm
res2 <-cv.glmnet(x, y, family="binomial", trace=TRUE)
my_lambda <- res2$lambda.1se
fit <- glmnet(x,y, family="binomial", lambda=my_lambda)

relaxed_res2 <- relax.glmnet(fit, x=x, y=y, 
                             family="binomial", path=TRUE)
par(mfrow=c(1,3))
plot(relaxed_res2, gamma=1)# default
plot(relaxed_res2, gamma=0.5)
plot(relaxed_res2, gamma=0)
relaxed_res2 <- relax.glmnet(fit, x=x, y=y, family="binomial", path=TRUE)
relaxed_res_gamma0 <- relax.glmnet(fit, x=x, y=y, gamma=0, family="binomial", path=TRUE)
relaxed_res_gamma0.5 <- relax.glmnet(fit, x=x, y=y, gamma=0.5, family="binomial", path=TRUE)
relaxed_res_gamma1 <- relax.glmnet(fit, x=x, y=y, gamma=1, family="binomial", path=TRUE)
# why aren't these the same? gamma = 1 should be lambda

relaxed_res_gamma0 <- relax.glmnet(fit, x=x, y=y, gamma=c(0, 0.25, 0.5, 0.75, 1), family="binomial", path=TRUE)
relaxed_res_gamma0.5 <- relax.glmnet(fit, x=x, y=y, gamma=0.5, family="binomial", path=TRUE)


pred <- predict(res, newx=x)
assess.glmnet(res, newx=x, newy=y) # results using lambda.1se
# grab some sex labeling data and look at it



# if (relax){
#   fit <- glmnet(X_train2[,rand_genes],Y_train2, family="binomial", lambda=my.lambda, alpha=my.alpha)
#   relaxo <- relax.glmnet(fit, x=X_train2[,rand_genes], y=Y_train2, gamma=0, family="binomial", path=TRUE, trace=TRUE)
#   relax0.5 <- relax.glmnet(fit, x=X_train2[,rand_genes], y=Y_train2, gamma=0.5, family="binomial", path=TRUE, trace=TRUE)
#   relax1 <- relax.glmnet(fit, x=X_train2[,rand_genes], y=Y_train2, gamma=1, family="binomial", path=TRUE, trace=TRUE)
#   
#   mat_coef <- coef(relaxo) %>% as.matrix()
#   nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
#   coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef))
#                              
#   print(nrow(coef_df))
#   mat_coef <- coef(relaxo) %>% as.matrix()
#   nonzero_coef2 <- data.frame(mat_coef[mat_coef[,1]!=0,])
#   print(nrow(nonzero_coef2))
#   
# }
# do we get better performance if we relax our results??
# IDK.
