
require('tidyverse')
metadata <- read.csv(sprintf("data/%s_metadata2.csv", prefix))
no_cell <- metadata %>% filter(cl_line=="")
refine_mapped <- read_csv(sprintf("data/%s_acc_to_cl.csv", prefix))
cell <- refine_mapped 
# // TODO - need to dedup cell_acc
both_cell <- no_cell %>% anti_join(cell, by=c("acc"="gsm")) %>% 
  select(acc) %>%  mutate(cl="", accession="") %>%
  rename(gsm=acc) %>% select(accession, cl, gsm) %>%
  bind_rows(cell) %>%
  rename(cl_descript_sample=cl, cl_acc_sample=accession) %>%
  mutate(is_sample_cl=ifelse(cl_acc_sample=="" & cl_descript_sample=="", FALSE, TRUE))

# we also want the metadata to -NOT- have CL mentions
cl_stud <- fread("../labeling/geo2drug/data/02_labeled_data/cell_line_mapped_gse.txt")

mapping <- fread("data/human_exp_to_sample.csv", data.table = FALSE)
study_sample_cl <- inner_join(cl_stud, mapping, by=c("gse"="study_acc"))
study_cl2 <- study_sample_cl %>% select(gse, sample_acc, everything()) %>% 
  rename(cl_acc_study=accession, cl_descript_study=cl, 
         cl_mention_study=cell_line, gsm=sample_acc) %>%
  mutate(is_study_cl=ifelse(cl_acc_study=="NA" & 
                              cl_descript_study=="NA" & 
                              !cl_mention_study, FALSE, TRUE))

# randomly sample some of these

cl_annot_all <- inner_join(both_cell, study_cl2, by="gsm") %>% 
  inner_join(human_metadata %>% select(acc, f_idx, idx), by=c("gsm"="acc")) %>%
  select(gsm, is_sample_cl, is_study_cl, f_idx, idx) 

both_neg <- cl_annot_all %>% filter(!is_sample_cl & !is_study_cl) %>% unique()
both_pos <- cl_annot_all %>% filter(is_sample_cl & is_study_cl) %>% unique()

set.seed(304)
neg_ex <- sample(1:nrow(both_neg), 6000, replace=FALSE)
pos_ex <- sample(1:nrow(both_pos), 6000, replace=FALSE)
not_cell_train <- both_neg[neg_ex,]
cell_train <- both_pos[pos_ex,]



pos_samp <- cell_train %>% rename(acc=gsm) %>% select(acc, f_idx, idx) %>% 
  arrange(idx) %>% unique() 
pos_samp %>% write_csv("cell_pos_sample.csv")
neg_samp <- not_cell_train %>% rename(acc=gsm) %>% select(acc, f_idx, idx) %>% 
  arrange(idx) %>% unique() 
neg_samp %>% write_csv("cell_neg_sample.csv")
