# add count data + present info for RNA-seq

library(tidyverse)

list_samples <- read_csv("data/01_sample_lists/list_samples.csv") 
list_studies <- read_csv("data/01_sample_lists/list_studies.csv") 
samp_to_study <- read_csv("data/01_sample_lists/sample_to_study.csv")


rc <- read_csv("data/02_metadata/qc/all_human_read_counts.csv") %>%
  bind_rows(read_csv("data/02_metadata/qc/all_mouse_read_counts.csv"))
present <- read_csv("data/02_metadata/qc/human_samples_present.csv") %>%
  bind_rows(read_csv("data/02_metadata/qc/mouse_samples_present.csv"))

colnames(present) <- "sample_acc"
present2 <- present %>% mutate(present=TRUE)

list_samples2 <- list_samples %>% 
  left_join(rc, by="sample_acc") %>%
  left_join(present2, by="sample_acc") 

nrow(list_samples2)
summary(list_samples2 %>% select(read_count, present, scrna)) # appears alll present. HM

studies <- read_csv("data/02_metadata/qc/human_studies_missing_samples.csv") %>%
  bind_rows(read_csv("data/02_metadata/qc/mouse_studies_missing_samples.csv"))
colnames(studies) <- "study_acc"



list_studies %>% semi_join(studies)  # 167
samp_to_study %>% semi_join(studies) # 762

# write this out

list_samples2 %>% write_csv("data/01_sample_lists/list_samples_w_qc.csv")

