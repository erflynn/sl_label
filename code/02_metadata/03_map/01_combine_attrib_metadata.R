# Code for combining the attribute metadata
# 
# In the process, gets counts that are used for descriptions in the paper.


require('tidyverse')
source("code/01_metadata/03_map/00_mapping_utils.R")


# --- load all the data --- #
ae_attrib <- read_csv("data/ae_attributes.csv")
gsm_attrib <- read_csv("data/gsm_key_value.csv")

load("data/sample_to_attr_sm.RData") # 137 MB --> sample_dat3
sample_dat4 <- sample_dat3 %>% 
  filter(key!="biomaterial provider")

# --- put it all together --- #

run_to_sample <- read_csv("data/sra_run_to_sample.csv")
rnaseq_dat <- sample_dat4 %>% 
  inner_join(run_to_sample, by=c("sample"="samples")) %>%
  rename(sample_acc=run) %>%
  select(sample_acc, key, value)

all_attrib <- rnaseq_dat %>%
  bind_rows(gsm_attrib %>% rename(sample_acc=gsm)) %>%
  bind_rows(ae_attrib) 

# NOTE - slow
all_attrib_clean <- all_attrib %>%
  mutate(across(c(key, value), clean_str, .names="{col}_clean")) 
all_attrib_clean <- all_attrib_clean %>% select(-study_acc)

# --- write it out --- #
save(all_attrib_clean, file="data/01_metadata/all_sample_attrib_clean.RData")


# ----- check counts! ----- #


# how many samples have *ANY* data
ae_files <- comb_metadata %>% 
  filter(!str_detect(sample_acc, "GSM|SRR|ERR|DRR") ) 
length(unique(ae_files$sample_acc))
#length(setdiff(unique(ae_files$sample_acc), unique(ae_attrib$sample_acc)))
# all AE files - note we're grabbing extra

microarray_s <- comb_metadata %>% 
  filter(str_detect(sample_acc, "GSM"))
length(microarray_s$sample_acc) # 448827
length(unique(gsm_attrib$gsm)) 
length(setdiff(microarray_s$sample_acc, unique(gsm_attrib$gsm))) 
# 5216 missing
# 448827 vs 453787

rnaseq_s <- comb_metadata %>% 
  filter(str_detect(sample_acc, "SRR|ERR|DRR"))
length(rnaseq_s$sample_acc) # 588,931
length(unique(rnaseq_dat$sample_acc))  # 568,430
length(setdiff(rnaseq_s$sample_acc,unique(rnaseq_dat$sample_acc)))
# 20502 hmm - are these the ones missing and included in MetaSRA?
head(setdiff(rnaseq_s$sample_acc,unique(rnaseq_dat$sample_acc)))

# how many of these are in comb_metadata and kept?
kept_s <- rnaseq_s %>% 
  filter(present & num_reads > 100000) # 241691
length(setdiff(kept_s$sample_acc,unique(rnaseq_dat$sample_acc))) # 9232
length(intersect(kept_s$sample_acc,unique(rnaseq_dat$sample_acc))) # 232459

head(setdiff(kept_s$sample_acc, unique(rnaseq_dat$sample_acc)))
# why are these missing? are they missing in the metadata?
# XX fraction are missing
# YY fraction overlap w the MetaSRA samples that are missing in this data?

