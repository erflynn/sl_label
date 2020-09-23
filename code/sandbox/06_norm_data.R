require('limma')
require('tidyverse')

load("data/pos_neg_genes.RData")
head(pos2[,1:5])
head(neg2[,1:5])

df <- cbind(pos2, neg2)
df

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", col_types="cccccccdld")
rc <- comb_metadata %>% 
  filter(data_type=="rnaseq", organism=="human") %>%
  select(sample_acc, num_reads)

get_count_mat <- function(dat, rc){
  dat_t <- data.frame(t(dat %>% select(-rid)))
  colnames(dat_t) <- dat$rid
  dat_t$sample_acc <- rownames(dat_t)
  dat_t2 <- dat_t %>% left_join(rc, by="sample_acc")
  count_matrix <- dat_t2 %>%
    mutate(across(dat$rid, ~.*num_reads/100000))
  count_df <- data.frame(t(count_matrix %>% select(-sample_acc, -num_reads))) 
  colnames(count_df) <- dat_t$sample_acc
  rownames(count_df) <- dat$rid  
  return(count_df)
}

pos_count_mat <- get_count_mat(pos2, rc) 
neg_count_mat <- get_count_mat(neg2, rc)

# look at the data...

# NOW TRY LIMMA VOOM
require('limma')
require('edgeR')

pos_count_mat2 <- pos_count_mat[,colSums(pos_count_mat) > 1000000]
neg_count_mat2 <- neg_count_mat[,colSums(neg_count_mat)> 1000000]
sex_lab <- c(rep(1, ncol(pos_count_mat2)), 
             rep(2, ncol(neg_count_mat2)))
design <- data.frame(sex_lab)
counts <- cbind(pos_count_mat2, neg_count_mat2)
dge <- DGEList(counts=counts)
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, plot=TRUE)
v2 <- voom(counts, design, plot=TRUE)

fit <- lmFit(v, data.frame(design))
fit <- eBayes(fit)
tt <- topTable(fit, coef=ncol(design), n=nrow(v$E)) 
tt$transcript <- rownames(tt)
tt2 <- tt %>% as_tibble() %>% left_join(xy_rnaseq) %>% 
  filter(adj.P.Val < 0.01 & abs(logFC) > 0.3)
tt2 %>% filter(logFC > 0) %>% arrange(adj.P.Val) 
tt2 %>% filter(logFC < 0) %>% arrange(adj.P.Val) 

pos_genes <- tt2 %>% filter(logFC > 0) %>% pull(transcript)
neg_genes <- tt2 %>% filter(logFC < 0) %>% pull(transcript)


# plot a couple things

# sample some pts
# see how it looks
num_f <- ncol(pos_count_mat2)
num_m <- ncol(neg_count_mat2)
rand_idx <- c(sample(1:num_f, 20), sample(1:num_m, 20)+num_f)
initialE <- counts[keep,rand_idx]
smallE <- v$E[,rand_idx]
boxplot(smallE)

v_norm <- normalize.quantiles(v$E)
rownames(v_norm) <- rownames(v$E)
posSums <- colSums(v_norm[pos_genes,])
negSums <- colSums(v_norm[neg_genes,])
sumG <- tibble("sex"=sex_lab,
  "pos"=posSums,
  "neg"=negSums)
normalize.quantiles()
ggplot(sumG, aes(x=pos, y=neg, col=factor(sex)))+geom_point(alpha=0.5)



pcs1 <- prcomp(t(smallE[pos_genes,]), center=FALSE)
pcs2 <- prcomp(t(smallE[neg_genes,]), center=FALSE)

plot(pcs1$x[,1], pcs2$x[,1], col=c(rep("blue", 20), rep("red",20)))
# we don't see anything

# grab X and Y chrosome genes

xy_rnaseq <- read_csv("data/human_xy_rnaseq.csv")
table(xy_rnaseq$chromosome_name)
xy_2 <- xy_rnaseq %>% filter(transcript %in% rownames(v$E))
xy_2 %>% filter(chromosome_name=="X") # 2,359
xy_2 %>% filter(chromosome_name=="Y") # 21, we lose a lot of Ychr genes

# ----- TRY TRAIN/TEST SPLIT ----- #

load("data/07_model_dat/human_rnaseq_sex_train_mat.RData")
train_data <- rownames(X_train)
test_data <- rownames(X_test)

pos_train <- pos_count_mat2[,intersect(colnames(pos_count_mat2), train_data)]
neg_train <- neg_count_mat2[,intersect(colnames(neg_count_mat2), train_data)]
pos_test <- pos_count_mat2[,intersect(colnames(pos_count_mat2), test_data)]
neg_test <- neg_count_mat2[,intersect(colnames(neg_count_mat2), test_data)]

sex_lab_train <- c(rep(0, ncol(pos_train)), 
             rep(1, ncol(neg_train)))
sex_lab_test <- c(rep(0, ncol(pos_test)), 
             rep(1, ncol(neg_test)))
design_train <- data.frame(sex_lab_train)
counts_train <- cbind(pos_train, neg_train)
dge_train <- DGEList(counts=counts_train)
keep_train <- filterByExpr(dge_train, design_train)
dge_train <- dge_train[keep_train,,keep.lib.sizes=FALSE]
dge_train <- calcNormFactors(dge_train)

v_train <- voom(dge_train, plot=TRUE) # normalize="quantile"
train_norm <- v_train$E
rowMeans <- apply(v_train$weights, 1, mean)

# normalize the test data the same way
counts_test <- cbind(pos_test, neg_test)
dge_test <- DGEList(counts=counts_test[keep_train,])
dge_test <- calcNormFactors(dge_test)
v_test <- voom(dge_test, plot=TRUE)
test_norm <- v_test$E

#devtools::install_github("LXQin/precision") #
require('PRECISION')



counts <- counts_train[keep_train,]
lib.size <- colSums(counts)
y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
train_norm2 <- y
qn <- quant.norm(train=train_norm2, test=NULL)
train_norm3 <- qn$train.qn

counts <- counts_test[keep_train,]
lib.size <- colSums(counts)
y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
test_norm2 <- y
test_qn <- quant.norm(test=test_norm2, ref.dis=qn$ref.dis)
test_norm3 <- test_qn$test.fqn

# how do we make sure this is ok for single sex studies?
#   samples in ss studies have the same distribution, right?


names(sex_lab_train) <- colnames(train_norm3)
sel_genes <- intersect(c(pos_genes, neg_genes), rownames(train_norm3))
cv.fit <- cv.glmnet(x=t(train_norm3[sel_genes,]), y=sex_lab_train, family="binomial", standardize=FALSE)
plot(cv.fit)
table(coef(cv.fit, s = "lambda.1se")[,1]!=0) # 113
train_pred <- predict(cv.fit, newx=t(train_norm3[sel_genes,]), s="lambda.1se", type="response", standardize=FALSE)
train_preds <- cbind("sample"=rownames(train_pred), "pred"=train_pred[,1],
                     "true_sex"=sex_lab_train) %>% as_tibble() %>%
  mutate(across(pred:true_sex, as.numeric))
train_preds2 <-train_preds %>% mutate(pred_sex=ifelse(pred > 0.5, 1, 0)) 
high_df <- train_preds2 %>% filter(true_sex!=pred_sex, pred > 0.90 | pred < 1) # 46

# remove these

sum(train_preds2$pred_sex==train_preds2$true_sex)/nrow(train_preds2) # 0.918
train_preds2 %>% head()


pred <- predict(cv.fit, newx=t(test_norm3[sel_genes, ]), s="lambda.1se", type="response", standardize=FALSE)
preds <- cbind("sample"=rownames(pred), "pred"=pred[,1], "true_sex"=sex_lab_test) %>% as_tibble()
preds2 <-preds %>% mutate(pred=ifelse(pred > 0.5, 1, 0))
sum(preds2$pred==preds2$true_sex)/nrow(preds2) # 0.9
