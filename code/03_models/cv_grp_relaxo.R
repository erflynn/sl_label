require('data.table')
require('tidyverse')
require('groupdata2')
require('glmnet')
#install.packages("~/Downloads/glmnet_3.0-2.tgz", repos=NULL)
# **NEED GLMNET3.0**

source("code/03_models/cv_utils.R")

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
data_type <- args[2]
ds <- args[3]

# set up the data
pos <- fread(sprintf("data/%s/%s/04_sl_input/%s_%s_pos_expr.csv", data_type, prefix, prefix, ds ), data.table=FALSE)
neg <- fread(sprintf("data/%s/%s/04_sl_input/%s_%s_neg_expr.csv", data_type, prefix, prefix, ds ), data.table=FALSE)

if (data_type == "compendia"){
  #mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample.csv", prefix), data.table=FALSE)  
  mapping <- fread(sprintf("data/01_metadata/%s_%s_dedup_mapping.csv", prefix, data_type), data.table = FALSE) %>%
    rename(study_acc=study_str)
} else {
  MIN.READS <- 100000
  # RNA-seq data does not have duplicates... interesting
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix), data.table=FALSE)  %>%
    filter(present & num_reads >= MIN.READS)
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
train_dat <- train_test_split %>% filter(partition!=5)
test_dat <- train_test_split %>% filter(partition==5)

# grab relevant info
fold_list <- partition_group_data(train_dat %>% select(-partition), nfolds=8)

# --- preprocess as needed (e.g. boxcox for RNASeq) --- #
if (data_type=="rnaseq"){
  
}

# CV folds


# for each fold in folds
#  for each param in params:
#    cv.glmnet() on all but fold using params (alpha)
#    predict on fold
#    return(pred_accuracy)