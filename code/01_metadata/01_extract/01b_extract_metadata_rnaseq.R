# grab the RNA-seq expression metadata and extract relevant terms
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
metadata_list <- fromJSON(file = sprintf("data/rnaseq/%s/aggregated_metadata.json", prefix))

# ----- 1. experiment level ----- #
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
  write.csv(file=sprintf("data/01_metadata/%s_rnaseq_experiment_metadata.csv", prefix), 
            quote=TRUE, row.names=FALSE)



# ----- 2. study to sample mapping ----- #

mapping <- exp.df2 %>% 
  as.data.frame() %>% 
  select(study_acc, sample_acc) %>%
  separate_rows(sample_acc, sep=";") %>%
  arrange(study_acc, sample_acc)
mapping %>% 
  write.csv(file=sprintf("data/01_metadata/%s_rnaseq_exp_to_sample.csv", prefix), 
            quote=TRUE, row.names=FALSE)


# ----- 3. sample level ----- #
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
sex_lab <- feat.df3 %>% filter(sex %in% c("male", "female")) 
cl.counts <- table(feat.df3$cl_line)
trt.counts <- table(feat.df3$trt)

feat.df3 %>% 
  write.csv(file=sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", prefix), 
            quote=TRUE, row.names=FALSE)


# // TODO arrange in some way with the mapping for indexes?

# grab the sex labeled data
sex_lab2 <- sex_lab %>% 
  select(acc, sex) %>% 
  left_join(mapping, by=c("acc"="sample_acc")) %>%
  rename(sample_acc=acc) %>%
  select(study_acc, sample_acc, sex) %>%
  arrange(study_acc, sample_acc)

sex_lab2 %>% write_csv(sprintf("data/%s/02_sample_lists/%s_rnaseq_sex_lab.csv", prefix,prefix))
