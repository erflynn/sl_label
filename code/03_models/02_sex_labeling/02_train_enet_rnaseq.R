
require('tidyverse')
require('glmnet')
require('data.table')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
idx <- as.numeric(args[2])

sl <- read_csv(sprintf("data/01_metadata/%s_rnaseq_sex_lab.csv", 
                       prefix, prefix))
exp_samp <- read_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample2.csv", prefix))
sl_present <- sl %>% semi_join(exp_samp %>% filter(present) %>% select(sample_acc))

miceadds::load.RData(sprintf("data/%s/04_sl_input/%s_rnaseq_sl_in.RData", prefix, prefix))
# then filter for what *is* present

# if they're not all the same length
#rownames(sample_expr_df) <- dat1$gene_name
rownames(sample_expr_df) <- sample_expr_df$gene_name


# load the transcript to gene mapping
map_l <- fread(sprintf("data/rnaseq/transcriptome_indices/%s/long/genes_to_transcripts.txt", prefix), 
               header=FALSE, data.table = FALSE)
#map_s <- fread(sprintf("data/rnaseq/transcriptome_indices/%s/short/genes_to_transcripts.txt", prefix),
#               header=FALSE, data.table=FALSE)
## map_s is the same as map_l
colnames(map_l) <- c("gene", "transcript")


# add rownames
# make sure the data are ok
load(sprintf("../sex_labeling/geo_pipeline/gpl_ref/%s/%s_gene_map.RData", prefix, prefix))
xy_genes <- ref_dat %>% 
  filter(chromosome_name %in% c("X", "Y")) %>% 
  select(chromosome_name, ensembl_gene_id) %>% 
  unique()

xy_transcripts <- inner_join(xy_genes, map_l, by=c("ensembl_gene_id"="gene"))

gene_names <- rownames(sample_expr_df)
overlap.transcripts <- xy_transcripts %>% filter(transcript %in% gene_names)
write_csv(overlap.transcripts, sprintf("data/%s/04_sl_input/xy_genes_rnaseq.csv", prefix))

sample_expr_df2 <- sample_expr_df[overlap.transcripts$transcript,] %>% select(-gene_name)
load("human_xy_expr_df.RData")

sample_expr_df2 <- sample_expr_df2[-which(duplicated(sample_expr_df2 )),]

# separate + shuffle
sl_sm <- sl %>% filter(sample_acc %in% colnames(sample_expr_df2))
sl <- data.frame(sl)
rownames(sl) <- sl$sample_acc
sl_sm <- sl[colnames(sample_expr_df2),]

# ------- do some plotting ---------- #
require('Rtsne')
dup_idx <- which(duplicated(data.frame(t(sample_expr_df2))))
cols2 <- colnames(sample_expr_df2)[-dup_idx]
sl_sm <- sl[cols2,]

# # there is a duplicated row -AND- column
# df3 <- data.frame(t(sample_expr_df2)) %>% distinct()
# tsne <- Rtsne(df3, dims = 2, 
#               perplexity=100, 
#               verbose=TRUE, max_iter = 1000)
# df <- data.frame(cbind(tsne$Y, sl_sm$sex))
# colnames(df) <- c("x1", "x2", "sex")
# df2 <- df %>% mutate(x1=as.numeric(x1), x2=as.numeric(x2))
# ggplot(df2, aes(x=x1, y=x2))+geom_point(aes(col=factor(sex)))
# 
# # look at this, this is confusion
# # // TODO - does the expression data look just as bad?
# pcs <- prcomp(sample_expr_df2, scale.=TRUE, center=TRUE)
# pc_df <- data.frame(cbind(pcs$rotation[,1:2], sl[colnames(sample_expr_df2),]$sex))
# colnames(pc_df) <- c("pc1", "pc2", "sex")
# pc_df2 <- pc_df %>% mutate(pc1=as.numeric(pc1), pc2=as.numeric(pc2))
# ggplot(pc_df2, aes(x=pc1, y=pc2))+geom_point(aes(col=factor(sex)))

# -------- look at the identity of things ------- #
sample_expr_df3 <- sample_expr_df2[,-dup_idx]
exp_sample <- read.csv("data/01_metadata/human_rnaseq_exp_to_sample2.csv")
metadata <- read.csv("data/01_metadata/human_rnaseq_sample_metadata.csv")
metadata_sm <- metadata %>% filter(acc %in% colnames(sample_expr_df3))
no_cl <- metadata_sm %>% filter(cl_line=="")
keep.cols <- colnames(sample_expr_df3) %in% no_cl$acc
sample_expr_df4 <- sample_expr_df3[,keep.cols]
sl_sm2 <- sl_sm[keep.cols,]

# make sure we divide such that each study is in a different fold..
study_counts <- sl_sm2 %>% group_by(study_acc) %>% count() %>% arrange(desc(n))
study_counts$group <- c(rep(c(1,2,3,1,2), nrow(study_counts)%/%5), rep(3, nrow(study_counts)%%5))

# foldid - vector of values between 1 and nfold 
# identifying which fold the obs is in

train_valid_data <- 
  study_counts %>% filter(group %in% c(1,3)) %>% 
  ungroup() 
train_valid_data2 <- train_valid_data %>%
  mutate(group2=c(rep(c(1:14), nrow(train_valid_data)%/%14), c(1:(nrow(train_valid_data)%%14))))

test_data <- study_counts %>% filter(group==2)

# -------- train/test split --------- #

# -------- CV split?? ------ #

outer_cv <- 
inner_cv <- 

sample_expr_df2 <- sample_expr_df4
set.seed(310)
f_df_u <- sample_expr_df2[,sl_sm2$sex=="female"]
f_df <- f_df_u[,sample(1:ncol(f_df_u), ncol(f_df_u)) ]
m_df_u <- sample_expr_df2[,sl_sm2$sex=="male"]
m_df <-m_df_u[,sample(1:ncol(m_df_u), ncol(m_df_u)) ]
  
# train/test
expr_f_train <- f_df[,1:ceiling(ncol(f_df)/2)]
expr_f_test <- f_df[,(ceiling(ncol(f_df)/2)+1):ncol(f_df)]

expr_m_train <- m_df[,1:ceiling(ncol(m_df)/2)]
expr_m_test <- m_df[,(ceiling(ncol(m_df)/2)+1):ncol(m_df)]

expr_train <- cbind(expr_f_train, expr_m_train)
train_sex_lab <- c(rep(0, ncol(expr_f_train)), rep(1, ncol(expr_m_train)))

expr_test <- cbind(expr_f_test, expr_m_test)
test_sex_lab <- c(rep(0, ncol(expr_f_test)), rep(1, ncol(expr_m_test)))

# train model
# train lasso on these data
x_train <- as.matrix(t(expr_train))
x_train <- apply(x_train, c(1,2),as.numeric)

x_test <- as.matrix(t(expr_test))
x_test <- apply(x_test, c(1,2),as.numeric)

#apply(x_train, 2, mean)
#apply(x_train, 2, sd)

#means <- apply(sample_expr_df %>% 
#                 select(-gene_name), 2, function(x) mean(x, na.rm=TRUE))
#sds <- apply(sample_expr_df %>% 
#                 select(-gene_name), 2, function(x) sd(x, na.rm=TRUE))

# all the means are the same - is this because of the TPM count?
# sds are not -- why do some have a higher standard deviation?


# ---------- denoising --------- #
# require("NoiseFiltersR")
# sample2 <- x_train %>% as_data_frame() %>% mutate(sex=train_sex_lab) 
# sample2 <- sample2 %>% mutate(sex=as.factor(sex))
# enn_res <- ENN(sex ~ ., data=sample2, k=5)
# 
# enn_res$remIdx

# remove most extreme...
cvfit = cv.glmnet(x_train, train_sex_lab, 
                  family="binomial", 
                  alpha=1, 
                  standardize=TRUE, 
                  type.measure="class")
preds_train <- predict(cvfit, newx=x_train, s="lambda.1se", type="response")
preds_class_train <- sapply(predict(cvfit, newx=x_train, s="lambda.1se", type="class"), as.numeric)
 sum(preds_class_train==train_sex_lab)/length(train_sex_lab) 

misl <- which(preds_class_train != train_sex_lab & (preds_train > 0.8 | preds_train < 0.2))

x_train2 <- x_train[-misl,]
train_sex_lab2 <- train_sex_lab[-misl]

cvfit2 = cv.glmnet(x_train2, train_sex_lab2, 
                  family="binomial", 
                  alpha=1, 
                  standardize=TRUE, 
                  type.measure="class")
preds_train <- predict(cvfit2, newx=x_train2, s="lambda.1se", type="response")
preds_class_train <- sapply(predict(cvfit2, newx=x_train2, s="lambda.1se", type="class"), as.numeric)

misl2 <- which(preds_class_train != train_sex_lab2 & (preds_train > 0.8 | preds_train < 0.2))

x_train3 <- x_train2[-misl2,]
train_sex_lab3 <- train_sex_lab2[-misl2]

# POSSIBLE WAYS TO IMPROVE
# 1. filter out cell line data
# 2. filter out mislabeled
# 3. make sure training/testing are different
# 3. select alpha
# 4. nested CV
# 5. prcomp / pc regression
# 6. prcomp/tsne viz
# 7. something for working with mislabeled training data?

# --------- train model ------- #
list.alphas <- seq(0,1,0.1)
hyperparam_res2 <- lapply(list.alphas, function(my.alpha){
  cvfit = cv.glmnet(x_train, train_sex_lab, 
                    family="binomial", 
                    alpha=my.alpha, 
                    standardize=TRUE)
  preds_train <- predict(cvfit, newx=x_train, s="lambda.1se", type="response")
  preds_class_train <- sapply(predict(cvfit, newx=x_train, s="lambda.1se", type="class"), as.numeric)
  train_acc <- sum(preds_class_train==train_sex_lab)/length(train_sex_lab) 
  
  preds_test <- predict(cvfit, newx=x_test, s="lambda.1se", type="response")
  preds_class_test <- sapply(predict(cvfit, newx=x_test, s="lambda.1se", type="class"), as.numeric)
  test_acc <- sum(preds_class_test==test_sex_lab)/length(test_sex_lab) 
  print(my.alpha)
  print(train_acc)
  print(test_acc)
  return(list("alpha"=my.alpha, "train_acc"=train_acc, "test_acc"=test_acc))
})

# which are likely mislabeled?


#[1] 0.8214423 0.9089845 0.8391188 0.1297520
names(preds_train)[preds_class_train != train_sex_lab & (preds_train < 0.6 & preds_train > 0.4)]

preds_test[preds_class_test != test_sex_lab & (preds_test > 0.8 | preds_test < 0.2)]
plot(density(preds_train[preds_class_train != train_sex_lab,]))

plot(density(preds_test[preds_class_test != test_sex_lab,]))
read_counts <- read_csv("data/human_read_counts_0.csv")
x_train2 <- x_train[intersect(read_counts$sample_acc, rownames(x_train)),]
x_test2 <- x_test[intersect(read_counts$sample_acc, rownames(x_test)),]
# boo only 23


mat_coef <- coef(cvfit, lambda="lambda.1se") %>% as.matrix()
nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef))
coef_df2 <- coef_df %>% filter(gene!="(Intercept)") %>%  mutate(coef=as.numeric(as.character(coef))) %>%
  arrange(coef) 
coef_df2 <- coef_df %>% left_join(xy_transcripts, by=c("gene"="transcript"))

in_dat <- x_train[,coef_df2$coef]
require('mclust')
mod2 <- MclustDA(in_dat, train_sex_lab)
summary(mod2)
plot(mod2, what = "scatterplot")

coef_df2$gene[which.min(coef_df2$coef)]
coef_df2$gene[which.max(coef_df2$coef)]

coef_df2 %>% 
  mutate(coef=as.numeric(as.character(coef))) %>%
  arrange(coef) %>% 
  filter(!is.na(chromosome_name)) %>% 
  write_csv(sprintf("data/%s/04_sl_input/rnaseq_%s_coef.csv", prefix, prefix))

save(cvfit, file=sprintf("data/%s/04_sl_input/rnaseq_cvfit.RData", prefix))

