
prefix <- "mouse"
data_type <- "rnaseq"
ds <- "sex"


require('tidyverse')
full_df <- do.call(rbind, lapply(3:6, function(idx){
  load(sprintf("data/06_fold_dat/fold_%s_%s_%s_%s.RData", prefix, data_type, ds, idx))
  df <- do.call(rbind, lapply(res, function(x) x$tv))
  df$fold <- idx
  return(df)
}))

full_df %>% write_csv(sprintf("data/06_fold_dat/%s_%s_%s_combined_params.csv", prefix, data_type, ds))



#full_df <- read_csv("data/mouse_microarray_sex_combined_params.csv") 
full_df2 <- full_df %>%
  rename(deviance="deviance.1",
         class="class.1",
         mse="mse.1",
         mae="mae.1")
full_long <- full_df2 %>% pivot_longer(cols=c("deviance","class","auc","mse","mae"), names_to="metric")

ggplot(full_long %>% filter(grp=="valid" & adaptive==FALSE), 
       aes(y=value, x=factor(alpha)))+
  geom_point(aes(col=metric), position=position_jitter(0.1))

ggplot(full_long %>% filter(metric=="class" & grp=="valid"), 
       aes(y=value, x=factor(alpha), fill=adaptive))+
  geom_boxplot()
ggplot(full_long %>% filter(metric=="class" & adaptive==FALSE), 
       aes(y=value, x=factor(alpha), fill=grp))+
  geom_boxplot()

ggplot(full_long %>% filter(metric=="deviance" & adaptive==FALSE), 
       aes(y=value, x=factor(alpha), fill=grp))+
  geom_boxplot()


ggplot(full_long %>% filter(metric=="class" & adaptive==FALSE), 
       aes(y=value, x=factor(alpha), fill=grp))+
  geom_boxplot()


ggplot(full_long %>% filter(metric=="auc" & adaptive==FALSE), 
       aes(y=value, x=factor(alpha), fill=grp))+
  geom_boxplot()



summarized <- full_long %>% filter(!adaptive &grp=="valid") %>%
  group_by(alpha, metric) %>%
  summarize(value_mu=median(value),
         value_sd=sd(value, na.rm=TRUE),
         lambda_mu=median(lambda),
         lambda_sd=sd(lambda, na.rm=TRUE))


 summarized %>% 
  filter(metric=="auc") %>% 
  ungroup() %>% 
  filter(value_mu==max(value_mu))

best_pars <- summarized %>% 
  filter(metric=="class") %>% 
  ungroup() %>% 
  filter(value_mu==min(value_mu))

best.alpha <- best_pars$alpha
best.lambda <- best_pars$lambda_mu

ggplot(full_df2 %>% filter(!adaptive & grp=="valid" & lambda < 0.5 ), 
       aes(x=lambda, y=class)) + geom_point(aes(col=factor(alpha)))


new_idx <- which(seq(0,1,0.1)==best.alpha)

all_ydf <- do.call(rbind, lapply(1:6, function(idx){
  load(sprintf("data/06_fold_dat/fold_%s_%s_%s_%s.RData", prefix, data_type, ds, idx))
  return(res[[new_idx]]$ydf)
}))

ydf2 <- all_ydf %>% mutate(ypred=round(as.numeric(yhat))) 

table(ydf2[,c("y","ypred")])

cv_accuracy <- sum(ydf2$y==ydf2$ypred)/nrow(all_ydf) 
# mouse microarray (all genes) 0.922 --> 0.920
# human microarray 0.912
# rat 0.814
# human rnaseq 0.876
# mouse rnaseq 0.935

######
require('data.table')
require('tidyverse')
require('groupdata2')
require('glmnet')
require('bestNormalize')
#install.packages("~/Downloads/glmnet_3.0-2.tgz", repos=NULL)
# **NEED GLMNET3.0**

source("code/03_models/cv_utils.R")

options(stringsAsFactors = FALSE)
NFOLDS <- 6

# set up the data

if (data_type == "microarray"){
  #mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample.csv", prefix), data.table=FALSE)  
  mapping <- fread(sprintf("data/01_metadata/%s_compendia_dedup_mapping.csv", prefix, data_type), data.table = FALSE) %>%
    rename(study_acc=study_str)
  pos <- fread(sprintf("data/05_train_df/%s_%s_%s_pos_expr.csv",prefix, data_type, ds ), data.table=FALSE)
  neg <- fread(sprintf("data/05_train_df/%s_%s_%s_neg_expr.csv", prefix, data_type, ds ), data.table=FALSE)
  
} else {
  MIN.READS <- 100000
  # RNA-seq data does not have duplicates... interesting
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix), data.table=FALSE)  %>%
    filter(present & num_reads >= MIN.READS)
  miceadds::load.Rdata(sprintf("data/05_train_df/%s_rnaseq_%s_pos_expr.RData", prefix, ds), "pos")
  miceadds::load.Rdata(sprintf("data/05_train_df/%s_rnaseq_%s_neg_expr.RData", prefix, ds), "neg")
  pos <- pos %>% rename(rid=gene_name)
  neg <- neg %>% rename(rid=gene_name)
}



pos_map <- mapping %>% filter(sample_acc %in% colnames(pos))
neg_map <- mapping %>% filter(sample_acc %in% colnames(neg))
my_dat <- bind_rows(pos_map %>% mutate(class=1), 
                    neg_map %>% mutate(class=0)) %>%
  mutate(class=factor(class), study_acc=factor(study_acc))

# training + testing division
set.seed(304)
train_test_split <- partition_group_data(my_dat, grp_col="study_acc", nfolds=5)
train_test_split %>% group_by(partition, class) %>% count()
train_part <- train_test_split %>% filter(partition!=5) %>% arrange(class, sample_acc)
test_part <- train_test_split %>% filter(partition==5) %>% arrange(class, sample_acc)

Y_train <- train_part$class
Y_test <- test_part$class

# divide into folds
fold_list <- partition_group_data(train_part %>% select(-partition), nfolds=NFOLDS) %>% arrange(class, sample_acc)
samp_to_fold <- fold_list$partition
names(samp_to_fold) <- fold_list$sample_acc

# grab relevant dfs
X_train <- t(cbind(neg[,train_part$sample_acc[Y_train==0]],pos[,train_part$sample_acc[Y_train==1]]))
X_test <- t(cbind(neg[,test_part$sample_acc[Y_test==0]],pos[,test_part$sample_acc[Y_test==1]]))
colnames(X_train) <- neg$rid
colnames(X_test) <- neg$rid
rownames(X_train) <- c(train_part$sample_acc)
rownames(X_test) <- c(test_part$sample_acc)

if (ds=="sex"){
  
  if (data_type=="rnaseq"){
    xy_genes <- read_csv(sprintf("data/%s/%s/03_model_in/xy_genes_%s.csv", data_type,prefix, data_type))
    xy_genes_l <- xy_genes$transcript
  } else {
    xy_genes <- read_csv(sprintf("data/%s/%s/04_sl_input/xy_genes.csv", data_type,prefix))
    xy_genes_l <- xy_genes$overlap.genes
  }
  X_test <- X_test[,xy_genes_l]
  X_train <- X_train[,xy_genes_l]
}
#####

if (data_type=="rnaseq"){
  X_train[is.na(X_train)] <- 0
  num_zeros <- apply(X_train, 2, function(x) sum(x==0)) # remove > 70% zeros
  X_train <- X_train[,(num_zeros/nrow(X_train) <= 0.7)]
  
  list_genes <- colnames(X_train)
  X_test <- X_test[,list_genes]
  
  X_train <- apply(X_train, 2, function(col) boxcox(col+0.5)$x.t) 
  X_test <- apply(X_test, 2, function(col) boxcox(col+0.5)$x.t)
  
}


fit <- glmnet(X_train, Y_train, family="binomial", 
              #lambda=best.lambda, 
              alpha=best.alpha,
              standardize=TRUE)


#preds_train <- predict(fit, newx=X_train,  type="response")
preds_class_train <- sapply(predict(fit, newx=X_train,  s=best.lambda, type="class"), as.numeric)
sum(preds_class_train==Y_train)/length(Y_train) 
# 0.930 --> 0.925 mouse microarray
# 0.916 (human microarray)
# 0.877 rat
# 0.878 (human RNA-seq)
# 0.935 (mouse RNA-seq)
preds_test <- predict(fit, newx=X_test,  s=best.lambda, type="response")
preds_class_test <- sapply(predict(fit, newx=X_test, s=best.lambda, type="class"), as.numeric)
sum(preds_class_test==Y_test)/length(Y_test) 
# 0.899 (human microarray)
# 0.896 --> 0.899 mouse micorarray
# 0.801 rat
# 0.882 (human RNA-seq)
# 0.913 (mouse RNA-seq)

my.lambda <- best.lambda
# human -r: 0.07086
# human -m - 0.03618 
# mouse -m 0.03909
mat_coef <- coef(fit, s=my.lambda) %>% as.matrix()
save(fit, my.lambda, file=sprintf("data/07_model_dat/fit_%s_%s_%s.RData", prefix, data_type, ds))
nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef))

coef_df2 <- coef_df %>% arrange(as.numeric(as.character(coef)))

load(sprintf("../sex_labeling/geo_pipeline/gpl_ref/%s_gene_map.RData", prefix))

if (prefix=="human"){
  ref_cols <- c("ensembl_gene_id", "hgnc_symbol",  "chromosome_name")
} else {
  ref_cols <- c("ensembl_gene_id",  "chromosome_name")
}

if (data_type=="rnaseq"){
  # how do we map these?
  coef_df3 <- coef_df2 %>% 
    rename(transcript=gene) %>%
    left_join(xy_genes %>% select(-chromosome_name), by=c("transcript")) %>%
    left_join(ref_dat %>% select(all_of(ref_cols)) %>% unique(), 
              by=c("ensembl_gene_id"))

} else {
  coef_df3 <- coef_df2 %>% 
    left_join(ref_dat %>% select(all_of(ref_cols)) %>% unique(), 
              by=c("gene"="ensembl_gene_id"))
}


coef_df3 %>% write_csv(sprintf("data/07_model_dat/%s_%s_%s_coef.csv", prefix, data_type, ds))
save(X_train, X_test, Y_train, Y_test, 
     file=sprintf("data/07_model_dat/%s_%s_%s_train_mat.RData", prefix, data_type, ds))

