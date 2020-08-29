# metadata_ss_hc.R
# E Flynn 
#
# Make lists of:
#  - single sex studies
#  - mixed sex studies 
# w/ >= 10 samples
#
# Setup code for HC testing:
#  - study-by-study: read in, cluster 
#  - try this on the new matrix

require('tidyverse')
MIN.READS <- 100000
prefix <- "human"
metadata_sex <- read.csv(sprintf("data/01_metadata/%s_rnaseq_metadata_sex.csv", prefix),
                         stringsAsFactors = FALSE) %>% as_tibble()

exp_sample <- read_csv(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix))

sex_lab <- metadata_sex %>% filter(sex %in% c("male", "female"))
sex_lab2 <- exp_sample %>% left_join(sex_lab %>% rename(sample_acc=acc))
study_samp2 <- sex_lab2 %>% mutate(qc_pass=(present & (num_reads >= MIN.READS))) %>%
  select(study_acc, sample_acc, mapped_sex, qc_pass) %>%
  rename(sex=mapped_sex)

study_counts_sex <- study_samp2 %>%
  group_by(study_acc) %>% 
  summarize( num_f=sum(sex=="female", na.rm=TRUE), 
             num_m=sum(sex=="male", na.rm=TRUE), num_samples=n(), 
             num_qc=sum(qc_pass)) %>%
  mutate(study_sex=case_when(
    num_samples < 10 ~ "small", # less than 10 samples
    num_f==num_samples ~ "female-only",
    num_m==num_samples ~ "male-only", 
    (num_m > 5 & num_f > 5) ~ "mixed sex",
    (num_f + num_m == 0) ~ "unknown",
    TRUE ~ "other")) # not all of the samples pass

table(study_counts_sex$study_sex)
mixed_sex_studies <- study_counts_sex %>% filter(study_sex=="mixed sex")
f_studies <- study_counts_sex %>% filter(study_sex=="female-only")
m_studies <- study_counts_sex %>% filter(study_sex=="male-only")

# do we have data on any of these that we can look at already?
load("~/Documents/EMILY/Stanford/Lab/AltmanLab/projects/drug_trt/data/human_xy_expr_df.RData")
# --> sample_expr_df2

# look at these
samples <- colnames(sample_expr_df2)
xy_genes <-read_csv("data/xy_genes_rnaseq.csv")
expr2 <- sample_expr_df2[xy_genes$transcript,]

mixed_sex_samples <- study_samp2 %>% semi_join(mixed_sex_studies) %>% unique()
f_samples <- study_samp2 %>% semi_join(f_studies) %>% unique()
m_samples <- study_samp2 %>% semi_join(m_studies) %>% unique()


# try clustering a few!

load("../labeling/gpl_ref/human_gene_map.RData")
xy_genes2 <- xy_genes %>% 
  left_join(ref_dat %>% select(ensembl_gene_id, hgnc_symbol) %>% unique())
toker_genes <- xy_genes2 %>% 
  filter(hgnc_symbol %in% c("XIST", "KDM5D", "RPS4Y1"))  
sl <- sex_lab %>% filter(acc %in% samples) %>% unique()


mm <- study_samp2 %>% filter(sample_acc %in% samples & sex %in% c("female", "male")) %>% unique()
mm_l <- mm %>% group_by(study_acc) %>% count() %>% filter(n>=10)
mm_filt <- mm %>% semi_join(mm_l)

require('limma')
require('rcompanion')
require('bestNormalize')
num_zeros <- apply(expr2, 1, function(x) sum(x==0)) # remove > 70% zeros

expr3 <- expr2[(num_zeros/ncol(expr2) <= 0.7),]
expr3 <- expr3[,unique(sl$acc)]
rowSums <- apply(expr3, 1, sum)
voom_trans <- voom(expr3) # assumes raw counts not TPM
blom_trans <- apply(expr3, 1, function(row) blom(row, alpha=3/8))
boxcox_trans <- apply(expr3, 1, function(row) boxcox(row+0.5)$x.t)
boxcox_obj <- apply(expr3[,sels], 1, function(row) boxcox(row+0.5))
boxcox_obj_t <- apply(expr3[,sels], 1, function(row) boxcox(row+0.5))
boxcox_train2 <- apply(expr3[,sels], 1, function(row) boxcox(row+0.5)$x.t)
boxcox_test2 <- apply(expr3[,-sels], 1, function(row) boxcox(row+0.5)$x.t)

expr3_test <- expr3[,-sels]

boxcox_res_test <- expr3_test
for (i in 1:nrow(expr3_test)){
  boxcox_res_test[i,] <- predict(boxcox_obj_t[[i]], newdata=(expr3_test[i,]+0.5))
} # SLOWWWW


expr4 <- t(blom_trans)
expr4 <- t(boxcox_trans)

expr4 <- voom_trans$E
# log scale everything! + normalize
# log_expr2 <- apply(expr2, c(1,2), function(x) log(x+1,2)) # slow
# norm_log_expr2 <- t(scale(t(log_expr2)))
# 
# mm_expr <- log_expr2[,unique(sl$acc)] # normalize each row?
toker_genes2 <- toker_genes %>% filter(transcript %in% rownames(expr4))
xist <- expr4[(toker_genes2 %>% filter(hgnc_symbol=="XIST"))$transcript,]
rps <- expr4[(toker_genes2 %>% filter(hgnc_symbol=="RPS4Y1"))$transcript,]
kdm <- expr4[(toker_genes2 %>% filter(hgnc_symbol=="KDM5D"))$transcript,]

df <- do.call(rbind, 
              list(apply(xist,2, sum), 
                   apply(rps,2, sum), 
                   apply(kdm,2, sum)))
df2 <- data.frame(cbind(t(df), "sex"=sl$mapped_sex)) 
colnames(df2) <- c("xist", "rps", "kdm", "sex")
df3 <- df2 %>% mutate(xist=as.numeric(as.character(xist)),
               rps=as.numeric(as.character(rps)),
               kdm=as.numeric(as.character(kdm)))

ggplot(df3, aes(x=xist, y=kdm))+geom_point(aes(col=sex))
ggplot(df3, aes(x=xist, y=rps))+geom_point(aes(col=sex))

# cluster and filter
res <- kmeans(df3[,1:3], 2)
pred_lab <- sapply(res$cluster, function(x) x-1)
true_lab <- sapply(df3$sex, function(x) ifelse(x=="female", 0, 1))
sum(pred_lab==true_lab)/length(true_lab) # 86.2%

mislab <- cbind(df3[which(pred_lab!=true_lab),], sample=mm_filt[which(pred_lab!=true_lab),]$sample_acc)
 # 96

all_z <- mislab %>% filter(rps < 0.1 & kdm < 0.1 & xist < 0.1) # 51

# which are actual flipped
flipped <- mislab %>% filter(
  (xist > 0.5 & rps < 0.5 & kdm < 0.5 & sex=="male") |
  ((xist < rps | xist < kdm) &  sex=="female")) # 8 potentially mislabeled or ambiguous

confusing <- mislab %>%  filter(!sample %in% all_z$sample & !sample %in% flipped$sample) # 37, all male
# appears some zeros here! :/
confusing %>% filter(rps > xist | kdm > xist) # 31/37 have one greater than xist

# train a model removing the flipped data
df4 <- df3 #[!(sl$acc %in% flipped$sample),]
sels <- sample(1:nrow(df4), 600)
mod <- glm(factor(sex) ~ xist + rps + kdm, data=df4[sels,], family="binomial")
pred_train <- predict(mod, newdata=df4[sels,], type="response")
pred_test <- predict(mod, newdata=df4[-sels,], type="response")
sum((sapply(pred_train, function(x) 
  ifelse(x < 0.5, "female", "male"))==df4[sels,]$sex))/length(sels) # 88.7%
sum((sapply(pred_test, function(x) 
  ifelse(x < 0.5, "female", "male"))==df4[-sels,]$sex))/(nrow(df4)-length(sels)) # 91.2%

# try lasso again with the data
require('glmnet')
sex_lab <- unlist(sapply(sl$sex, function(x) ifelse(x=="female", 0, 1)))
names(sex_lab) <- NULL
sels <- sample(1:nrow(df4), 1000)


# more filtering
head(mm_expr[,1:5])
num_zeros <- apply(mm_expr, 1, function(x) sum(x==0)) # remove > 70% zeros

mm_expr2 <- mm_expr[(num_zeros/ncol(mm_expr) <= 0.7),]
max_vals <- apply(mm_expr2, 1, max)
mean_vals <- apply(mm_expr2, 1, mean)
sd_vals <- apply(mm_expr2, 1, sd)
three_sd <-(mean_vals+3*sd_vals)
table(max_vals > six_sd)
mm_expr3 <- mm_expr2

# replace high values with the max for each row
for(i in 1:nrow(mm_expr3)){
  mm_expr3[i,][mm_expr3[i,] > three_sd[i]] <- three_sd[i]
}

max_vals3 <- apply(mm_expr3, 1, max)
table(max_vals3 > three_sd)

# filter flipped?

cvfit <- cv.glmnet(x=t(expr4[,sels]), y=factor(sex_lab[sels]) , family="binomial", standardize=FALSE)
preds_train <- predict(cvfit, newx=t(expr4[,-sels]), s="lambda.1se", type="response")
preds_class_train <- sapply(predict(cvfit, newx=t(expr4[,sels]), s="lambda.1se", type="class"), 
                            as.numeric)
preds_train[preds_class_train!=sex_lab[sels] & preds_train > 0.9]
preds_train[preds_class_train!=sex_lab[sels] & preds_train < 0.1]



sum(preds_class_train==sex_lab[sels])/length(sex_lab[sels])
preds_valid <- predict(cvfit, newx=t(expr4[,-sels]), s="lambda.1se", type="response")
preds_class_valid <- sapply(predict(cvfit, newx=t(expr4[,-sels]), s="lambda.1se", type="class"), as.numeric)
sum(preds_class_valid==sex_lab[-sels])/length(sex_lab[-sels]) 

# yay transformations!

cvfit <- cv.glmnet(x=boxcox_train2, y=factor(sex_lab[sels]) , family="binomial", standardize=FALSE)
preds_class_train <- sapply(predict(cvfit, newx=boxcox_train2, s="lambda.1se", type="class"), 
                            as.numeric)
sum(preds_class_train==sex_lab[sels])/length(sex_lab[sels])
preds_class_valid <- sapply(predict(cvfit, newx=boxcox_test2, s="lambda.1se", type="class"), as.numeric)
sum(preds_class_valid==sex_lab[-sels])/length(sex_lab[-sels]) 


# cv w transformed data?

mat_coef <- coef(cvfit, lambda="lambda.1se") %>% as.matrix()
nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef))
coef_df2 <- coef_df %>% left_join(xy_genes, by=c("gene"="transcript"))

coef_df2 %>%
  mutate(coef=as.numeric(as.character(coef))) %>%
  arrange(coef) %>% 
  filter(!is.na(chromosome_name)) 

# try running massiR

# ICA/PCA??

pc <- prcomp(expr4, scale=FALSE, center=FALSE)
pc_df <- data.frame(cbind(pc$rotation[,1:5], "sex"=sex_lab))
colnames(pc_df) <- c("pc1", "pc2","pc3","pc4","pc5", "sex")
pc_df2 <- pc_df %>% mutate(pc1=as.numeric(pc1), pc2=as.numeric(pc2))
ggplot(pc_df2, aes(x=pc1, y=pc2))+geom_point(aes(col=factor(sex)))


mod <- glm(sex ~ ., data=pc_df[sels,], family="binomial")
pred_train <- predict(mod, newdata=pc_df[sels,], type="response")
pred_test <- predict(mod, newdata=pc_df[-sels,], type="response")
sum((sapply(pred_train, function(x) 
  ifelse(x < 0.5, 0, 1))==pc_df[sels,]$sex))/length(pc_df[sels,]$sex) 
sum((sapply(pred_test, function(x) 
  ifelse(x < 0.5, 0,1))==pc_df[-sels,]$sex))/(length(pc_df[-sels,]$sex)) 

# USELESS


#### ----->
#  ANSWER:
#    blom or boxcox or voom
# then: glmnet
