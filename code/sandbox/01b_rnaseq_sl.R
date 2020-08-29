# generate the list of training data for sex labeling 
require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
MIN.READS <- 100000 

sl <- read_csv(sprintf("data/01_metadata/%s_rnaseq_sex_lab.csv", prefix, prefix))
mapping2 <- read_csv(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix))
mapping3 <- mapping2 %>% filter(present & num_reads >= MIN.READS)
sl2 <- sl %>% filter(sample_acc %in% mapping3$sample_acc)
sl2 %>% write_csv(sprintf("data/rnaseq/%s/03_model_in/%s_rnaseq_sl_filt.csv", prefix, prefix))

if (nrow(sl2) < 2000){
  exit()
}

# now sample to distribute across files as much as possible
set.seed(15)
sl3 <- sl2 %>% group_by(study_acc) %>% filter(n() >=5) %>% sample_n(5)
sl3.2 <- sl2 %>% group_by(study_acc) %>% filter(n() <5) 
sl4 <- rbind(sl3, sl3.2)

sl4 %>% ungroup() %>% 
  sample_n(2000) %>% 
  write_csv(sprintf("data/rnaseq/%s/03_model_in/%s_rnaseq_sl_training.csv", prefix, prefix))