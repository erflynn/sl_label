require('tidyverse')
require('data.table')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

# put together all the read counts into a file
list.f <- list.files(path=sprintf("data/rnaseq/%s/02_read_counts/", prefix))

count_df <- do.call(rbind, lapply(list.f, function(my.f)
  fread(sprintf("data/rnaseq/%s/02_read_counts/%s", prefix, my.f), data.table=FALSE)))

count_df %>% write_csv(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix))