# grab the microarray expression metadata and extract relevant terms
#
# this results in three files:
#   experiment_metadata
#   exp_to_sample
#   sample_metadata

library(rjson)
library(tidyverse)
require('data.table')

CHUNK.SIZE=20000

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

# ----- 0. grab metadata ------ #
metadata_list <- fromJSON(file = sprintf("data/%s/01_input/aggregated_metadata.json", prefix))
# experiment:
#   description, title, source_first_published, pubmed_id, technology, pubmed_id

list.feat <-
  lapply(metadata_list$samples, function(x)
    list("acc"=x$refinebio_accession_code,
         "sex"=x$refinebio_sex,
         "cl_line"=x$refinebio_cell_line,
         "compound"=x$refinebio_compound,
         "disease"=x$refinebio_disease,
         "age"=x$refinebio_age,
         "trt"=x$refinebio_treatment,
         "part"=x$refinebio_specimen_part,
         "src"=x$refinebio_source_database,
         "title"=x$refinebio_title,
         "platform"=x$refinebio_platform))

feat.df <- do.call(rbind, list.feat)
feat.df2 <- apply(feat.df , c(1,2), function(x) unlist(paste(x, sep=";")))
feat.df3 <- data.frame(feat.df2, stringsAsFactors = FALSE)

sex.counts <- table(feat.df3$sex)
# //TODO: this table is horrendous and needs to be reformatted
# --> unknown, mixed/pooled, female, male
sex_lab <- feat.df3 %>% filter(sex %in% c("male", "female")) # 67185 for mouse, 114856 for human
cl.counts <- table(feat.df3$cl_line)
trt.counts <- table(feat.df3$trt)

# sex_lab %>% 
#   filter(cl_line=="") %>% 
#   write_csv("sex_lab_no_cl.csv") # 45038

# ----- 1. how do sex_lab compare to the ALE sex lab? ------ #

ale_sex_lab <- read_csv(sprintf("../sex_labeling/geo_pipeline/data/01_sample_lists/%s_ale_sex_lab.csv", prefix))
feat.df3 <- feat.df3 %>% mutate(acc=as.character(acc), sex=as.character(sex))

length(intersect(ale_sex_lab$gsm, feat.df3$acc)) # 50,042 and 116,318 human
# mouse: there are 176,035 ale, and 279,998 in refinebio... hmm
# maybe this is because of differences in sample ids?
# -ALSO- all the ALE ones *HAVE* labels
ale2 <- ale_sex_lab %>% filter(gsm %in%  feat.df3$acc )
ale3 <- ale2 %>% left_join(feat.df3 %>% select(acc, sex), by=c("gsm"="acc"))
ale4 <- ale3 %>% mutate(ale_sex=case_when(
  Gender=="M" ~ "male", 
  Gender=="F" ~ "female",
  TRUE ~ Gender)) %>%
  rename(rb_sex=sex)
ale4 %>% filter(rb_sex==ale_sex & rb_sex!="") %>% nrow() 
# mouse: 37151, 74.2% agreement...
# human: 92094, 79.2% agreement

ale4 %>% filter(rb_sex!=ale_sex) %>% 
  select(rb_sex, ale_sex) %>% unique() 
# mouse: 78 diff mismatches 
# human: 438 diff mismatches

# BUT - there are only two exact swaps!
ale4 %>% filter(ale_sex=="male" & rb_sex=="female")
ale4 %>% filter(ale_sex=="female" & rb_sex=="male")
# GSE99192 - this is the study, and rb is wrong... (superseries + subseries)
# GSE17743, GSE55232, GSE83951 -- all of them rb says male and ale says female

# -----  create a training list ----- #

# label columns by file + filter for cols actually present
m_cols <- fread(sprintf("data/%s/02_sample_lists/%s_cols.txt",prefix, prefix), 
                data.table=FALSE, stringsAsFactors = FALSE)
m_to_idx <- data.frame(t(m_cols))
colnames(m_to_idx) <- c("col_name")
rownames(m_to_idx) <- NULL
m_to_idx <- m_to_idx %>% 
  mutate(idx=1:nrow(m_to_idx)) %>%
  mutate(f_idx=(idx-2)%/%CHUNK.SIZE) %>%
  filter(!is.na(col_name))

sex_lab2 <- sex_lab %>% inner_join(m_to_idx, by=c("acc"="col_name"))
sex_lab2 %>% 
  write_csv(sprintf("data/%s/02_sample_lists/%s_sex_lab_no_cl_idx.csv", prefix, prefix))

# randomly sample 1000 males and 1000 females
set.seed(228)
f_sample <- sex_lab2 %>% 
  filter(sex=="female") %>%
  sample_n(1000)
m_sample <- sex_lab2 %>% 
  filter(sex=="male") %>%
  sample_n(1000)
f_sample %>% write_csv(sprintf("data/%s/02_sample_lists/%s_f_sample.csv", prefix, prefix))
m_sample %>% write_csv(sprintf("data/%s/02_sample_lists/%s_m_sample.csv", prefix, prefix))


feat.df3 %>% 
  inner_join(m_to_idx, by=c("acc"="col_name")) %>%
  write.csv(file=sprintf("data/01_metadata/%s_metadata.csv", prefix), quote=TRUE, row.names=FALSE)


# ---------- grab study to sample mapping ---------- #
metadata_list <- fromJSON(file = sprintf("data/%s/01_input/aggregated_metadata.json", prefix))


list.exp <-
  lapply(metadata_list$experiments, function(x)
    list("study_acc"=unlist(x$accession_code),
         "title"=unlist(x$title),
         "description"=unlist(x$description),
         "date"=unlist(x$source_first_published),
         "sample_acc"=paste(x$sample_accession_codes, collapse=";")))

exp.df <- do.call(rbind, list.exp)
exp.df2 <- apply(exp.df , c(1,2), function(x) unlist(x))

exp_data <- exp.df2 %>% 
  as.data.frame() %>% 
  select(-sample_acc)
exp_data %>% 
  write.csv(file=sprintf("data/01_metadata/%s_experiment_metadata.csv", prefix), 
            quote=TRUE, row.names=FALSE)


# study to sample mapping 

mapping <- exp.df2 %>% 
  as.data.frame() %>% 
  select(study_acc, sample_acc) %>%
  separate_rows(sample_acc, sep=";") %>%
  arrange(study_acc, sample_acc)
mapping %>% 
  write.csv(file=sprintf("data/01_metadata/%s_exp_to_sample.csv", prefix), 
            quote=TRUE, row.names=FALSE)


