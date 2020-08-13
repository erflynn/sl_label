# DEFINE CELL LINE VS NON-CELL LINE, then set up train/test

require('tidyverse')
require('data.table')
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
data_type <- args[2]

if (data_type=="microarray"){
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample.csv", prefix, data_type), data.table = FALSE) 
  data_type2 <- "compendia"
  metadata <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", prefix))
} else {
  MIN.READS <- 100000
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix), data.table=FALSE)  %>%
    filter(present & num_reads >= MIN.READS)
  data_type2 <- data_type
  metadata <- read.csv(sprintf("data/01_metadata/%s_%s_sample_metadata.csv", prefix, data_type))
  
}


# drops samples such that no more than n samples are in a group
drop_greater_n <- function(ds, n_group, num_tot){
  ds2 <- ds %>% group_by(grp) %>% filter(n() >=n_group) %>% sample_n(n_group)
  ds3 <- ds %>% group_by(grp) %>% filter(n() < n_group) 
  ds4 <- rbind(ds2, ds3) %>% ungroup()
  if (nrow(ds4) <= num_tot){
    return(ds4)
  } else {
    ds5 <- ds4 %>% sample_n(num_tot)
    return(ds5)  
  }
}




sample_source_df <- read_csv(sprintf("data/02_labeled_data/%s_%s_sample_source.csv", prefix, data_type2))


set.seed(114)

pos_cl <- sample_source_df %>% filter(source_type=="named_cl") %>% left_join(mapping, by=c("acc"="sample_acc")) %>% select(-source_type)
neg_cl <- sample_source_df %>% filter(source_type=="tissue")  %>% left_join(mapping, by=c("acc"="sample_acc")) %>% select(-source_type)
# note - these do overlap!!


# randomly sample some of these
pos_cl_drop <- drop_greater_n(pos_cl %>% select(acc, study_acc) %>% dplyr::rename(grp=study_acc), 5, 2000) 
neg_cl_drop <- drop_greater_n(neg_cl %>% select(acc, study_acc) %>% dplyr::rename(grp=study_acc), 5, 2000)

if (data_type=="microarray"){
  pos_cl_drop <- pos_cl_drop %>% 
    left_join(metadata %>% select(acc, f_idx, idx)) %>% 
    select(-grp) %>%
    unique() %>%
    arrange(idx)
  neg_cl_drop <- neg_cl_drop %>% 
    left_join(metadata %>% select(acc, f_idx, idx)) %>% 
    select(-grp) %>%
    unique() %>%
    arrange(idx)
} else {
  pos_cl_drop <- pos_cl_drop %>% dplyr::rename(study_acc=grp, sample_acc=acc)
  neg_cl_drop <- neg_cl_drop %>% dplyr::rename(study_acc=grp, sample_acc=acc)
}

pos_cl_drop %>% sample_n(20) %>% left_join(mu_s3)
neg_cl_drop %>% sample_n(20) %>% left_join(mu_s3)

write_csv(neg_cl_drop, sprintf("data/03_train/%s_%s_cl_neg_sample.csv", prefix, data_type))
write_csv(pos_cl_drop, sprintf("data/03_train/%s_%s_cl_pos_sample.csv", prefix, data_type))
