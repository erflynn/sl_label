
require('tidyverse')
# get a list of training and testing data
load("data/07_model_dat/human_rnaseq_sex_train_mat.RData")
exp_to_sample <- read_csv("data/01_metadata/human_exp_to_sample_counts.csv")

train_df <- tibble("sample_acc"=rownames(X_train), "metadata_sex"=Y_train)
test_df <- tibble("sample_acc"=rownames(X_test), "metadata_sex"=Y_test)
train_df2 <- train_df %>% left_join(exp_to_sample, by="sample_acc") %>% 
  select(study_acc, everything())
test_df2 <- test_df %>% left_join(exp_to_sample, by="sample_acc") %>% 
  select(study_acc, everything())
# make sure it is one per study
set.seed(1)
train_df3 <- train_df2 %>% group_by(study_acc) %>% sample_n(1) %>% ungroup()
test_df3 <- test_df2 %>% group_by(study_acc) %>% sample_n(1) %>% ungroup()
# 573 train, 135 test

prefix <- "human"

# pull the count matrices --> full matrix of counts
sampleMat <- function(sample_id, study_id){
  sample.path <- sprintf("data/rnaseq/%s/00_infiles/%s/%s_quant.sf", 
                         prefix, study_id, sample_id);
  print(sample_id);
  if (file.exists(sample.path)){
    my_dat <- read_tsv(sample.path) 
    if (!"Name" %in% colnames(my_dat)){
      print(sprintf("corrupted file %s", sample.path))
      return(NA)
      #my_dat <- read.delim(sample.path, skip=1)
      #colnames(my_dat) <- c("Name", "Length","EffectiveLength","TPM", "NumReads")
    }
    quant_sf <- my_dat %>% dplyr::select(Name, TPM);
    colnames(quant_sf) <- c("gene_name", sample_id);
    
    #if (sample_id==sample_list$sample_acc[[1]]){
    #  return(quant_sf)
    #}
    
    return(quant_sf)
  } else {
    return(NA)
  }
}



# THEN try stuff
