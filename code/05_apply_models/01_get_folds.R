# look at the fold output from training models
require('tidyverse')
require('glmnet')
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
data_type <- args[2]
ds <- args[3]

# my_l <- list.files("data/06_fold_dat", pattern="fold_rat_microarray")
# 
# full_df <- do.call(rbind, lapply(my_l, function(my_f){
#   load(file=sprintf("data/06_fold_dat/%s", my_f))
#   df <- do.call(rbind, lapply(res, function(x) x$tv))
#   df$fold <- idx
#   return(df)
# }))

full_df <- do.call(rbind, lapply(1:6, function(idx){
  load(sprintf("data/06_fold_dat/fold_%s_%s_%s_%s.RData", prefix, data_type, ds, idx))
  df <- do.call(rbind, lapply(res, function(x) x$tv))
  df$fold <- idx
  return(df)
}))


# my_f <- sapply(1:8, function(it) {
#   lapply(1:4, function(idx){
#   my.f <- (sprintf("data/06_fold_dat/fold_%s_%s_%s_%s-it%s.RData", prefix, data_type, ds, idx, it))
#   if (file.exists(my.f)){
#     load(my.f)
#   } else {
#     return(NA)
#   }
#   df <- do.call(rbind, lapply(res, function(x) x$tv))
#   df$fold <- idx
#   df$it <- it
#   return(df)})
# })
# 
# my_df <- do.call(rbind, my_f[!is.na(my_f)])

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
  filter(metric=="class" & alpha!=0) %>% 
  ungroup() %>% 
  filter(value_mu==min(value_mu))

best.alpha <- best_pars$alpha
best.lambda <- best_pars$lambda_mu

ggplot(full_df2 %>% filter(!adaptive & grp=="valid" & lambda < 0.5 ), 
       aes(x=lambda, y=class)) + geom_point(aes(col=factor(alpha)))

best.ngenes <- 200 # //TODO
new_idx <- which(seq(0,1,0.1)==best.alpha)
ngenes_idx <- which(c(50, 100, 200, 300, 500)==best.ngenes)
all_ydf <- do.call(rbind, lapply(1:6, function(idx){
  load(sprintf("data/06_fold_dat/fold_%s_%s_%s_%s.RData", prefix, data_type, ds, idx))
  return(res[[new_idx]][[ngenes_idx]]$ydf)
}))

ydf2 <- all_ydf %>% mutate(ypred=round(as.numeric(as.character(yhat))) )

table(ydf2[,c("y","ypred")])

cv_accuracy <- sum(ydf2$y==ydf2$ypred)/nrow(all_ydf) 
# mouse microarray (all genes) 0.922 --> 0.920
# human microarray 0.912
# rat 0.814 m --> 0.737
# human rnaseq 0.876
# mouse rnaseq 0.935
# human microarray cl --> 0.813
# ---- use the fold output to train a new model ---- #

#save(X_train, X_test, Y_train, Y_test, 
load(file=sprintf("data/07_model_dat/%s_%s_%s_train_mat.RData", prefix, data_type, ds))

if (data_type=="microarray"){
  standardizeFlag = TRUE
} else {
  standardizeFlag = FALSE
}

if (best.ngenes < ncol(X_train)){
  require('limma')
  design <- as.matrix(as.numeric(as.character(Y_train)))
  fit <- lmFit(t(X_train), design)
  fit <- eBayes(fit)
  tt <-topTable(fit, number=ncol(X_train))
  tt$transcript <- rownames(tt)
  my_genes= head(tt$transcript, best.ngenes)
  X_train <- X_train[,my_genes]
  X_test <- X_test[,my_genes]
}

fit <- glmnet(X_train, Y_train, family="binomial", 
              #lambda=best.lambda, 
              alpha=best.alpha,
              standardize=standardizeFlag)


#preds_train <- predict(fit, newx=X_train,  type="response")
preds_class_train <- sapply(predict(fit, newx=X_train,  s=best.lambda, type="class"), as.numeric)
sum(preds_class_train==Y_train)/length(Y_train) 
# 0.930 --> 0.925 mouse microarray
# 0.916 (human microarray)
# 0.877 rat
# 0.878 (human RNA-seq)
# 0.935 (mouse RNA-seq)
# 0.777 rat microarray balanced
# 0.838

preds_test <- predict(fit, newx=X_test,  s=best.lambda, type="response")
preds_class_test <- sapply(predict(fit, newx=X_test, s=best.lambda, type="class"), as.numeric)
sum(preds_class_test==Y_test)/length(Y_test) 
# 0.899 (human microarray)
# 0.896 --> 0.899 mouse micorarray
# 0.801 rat
# 0.882 (human RNA-seq)
# 0.913 (mouse RNA-seq)
# 0.735 rat microarray balanced
# 0.813

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

