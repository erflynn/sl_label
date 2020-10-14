# code for aggregating count data into a single file
require('tidyverse')
args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
count_path = sprintf("data/rnaseq/%s/02_read_counts", organism)
count_f <- list.files(path=count_path, pattern="^count*")
count_dat <- map_df(count_f, function(x) 
  read_csv(sprintf("%s/%s", count_path, x), col_types="ccd")) %>% 
  bind_rows()
stopifnot(nrow(count_dat)==length(unique(count_dat$sample_acc)))

nrow(count_dat)
count_dat %>% filter(count > 100000) %>% nrow()

count_dat %>% write_csv(sprintf("%s/all_read_counts.csv", count_path))

