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

# // TODO: change file name to "compendia_sample" for consistency

feat.df3 %>% 
  inner_join(m_to_idx, by=c("acc"="col_name")) %>%
  write.csv(file=sprintf("data/01_metadata/%s_metadata.csv", prefix), quote=TRUE, row.names=FALSE)

# note - repeated
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

# // TODO: update the name to "compendia"
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


