require('glmnet')
require('tidyverse')
require('data.table')

prefix <- "human"

expr_cl <- fread(sprintf("data/%s/04_sl_input/cell_line_no_trt_expr.csv", prefix, prefix ), data.table=FALSE)
#expr_no_cl <- fread(sprintf("data/%s/04_sl_input/cell_neg_expr.csv", prefix, prefix ), data.table=FALSE)

cl_lab <- read_csv("data/cell_line_no_trt.csv")

rownames(expr_cl) <- expr_cl$rid
ord_expr_cl <- expr_cl[,cl_lab$acc]

set.seed(0317)
hela <- cl_lab %>% 
  filter(accession=="cvcl_0030") %>% 
  filter(acc %in% colnames(expr_cl)) %>%
  sample_n(nrow(hela)) # hela

non_hela <- cl_lab %>% filter(accession!="cvcl_0030")  %>% 
  filter(acc %in% colnames(expr_cl)) %>% sample_n(nrow(hela))

hela_expr <- expr_cl[,hela$acc]
non_hela_expr <- expr_cl[,non_hela$acc]
tot <- ncol(hela_expr)
midpt <- ceiling(tot/2)
hela_train <- hela_expr[,1:midpt]
hela_test <- hela_expr[,(midpt+1):tot]

non_hela_train <- non_hela_expr[,(midpt+1):tot]
non_hela_test <- non_hela_expr[,(midpt+1):tot]

expr_train <- cbind(hela_train, non_hela_train)
expr_test <- cbind(hela_test, non_hela_test)
train_lab <- c(rep(1, ncol(hela_train)), rep(0, ncol(non_hela_train)))
test_lab <- c(rep(1, ncol(hela_test)), rep(0, ncol(non_hela_test)))

x_train <- as.matrix(t(expr_train))
x_train <- apply(x_train, c(1,2),as.numeric)

x_test <- as.matrix(t(expr_test))
x_test <- apply(x_test, c(1,2),as.numeric)

# // set up nested CV
cvfit = cv.glmnet(x_train, train_lab, 
                  family="binomial", 
                  alpha=0.5, 
                  standardize=TRUE)
preds_train <- predict(cvfit, newx=x_train, s="lambda.1se", type="response")
preds_class_train <- sapply(predict(cvfit, newx=x_train, s="lambda.1se", type="class"), as.numeric)
sum(preds_class_train==train_lab)/length(train_lab) 

preds_test <- predict(cvfit, newx=x_test, s="lambda.1se", type="response")
preds_class_test <- sapply(predict(cvfit, newx=x_test, s="lambda.1se", type="class"), as.numeric)
sum(preds_class_test==test_lab)/length(test_lab) 

save(cvfit, file="data/hela_cvfit.RData")

plot(density(preds_test[preds_class_test != test_lab,]))
misl <- which(preds_class_test != test_lab & (preds_test > 0.6 | preds_test < 0.4))

# --- how does the choice of measure (e.g. "class", "dev" etc) reflect results?
#  need to understand this

# multinomial classifier
ord_expr_cl <- expr_cl[,cl_lab$acc]
ord_labels <- cl_lab$accession
set.seed(0339)
shuffled <- sample(1:length(ord_labels),ceiling(length(ord_labels)/2))

ord_expr_train <- ord_expr_cl[,shuffled]
ord_expr_test <- ord_expr_cl[,-shuffled]
ord_train_lab <- ord_labels[shuffled]
ord_test_lab <- ord_labels[-shuffled]
fit <- cv.glmnet(t(ord_expr_train), ord_train_lab, family="multinomial", alpha=0.5)
preds_class_train <- sapply(predict(cvfit, newx=t(ord_expr_train), s="lambda.1se", type="class"), as.numeric)
sum(preds_class_train==ord_train_lab)/length(ord_train_lab) 

preds_class_test <- sapply(predict(cvfit, newx=t(ord_expr_test), s="lambda.1se", type="class"), as.numeric)
sum(preds_class_test==ord_test_lab)/length(ord_test_lab) 

# what is the performance of a random classifier? ~ 0.10

