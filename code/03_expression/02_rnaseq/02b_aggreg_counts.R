# code for aggregating count data into a single file
require('tidyverse')
args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
count_path = sprintf("data/03_expression/rnaseq/%s/02_counts", organism)
count_f <- list.files(path=count_path, pattern="^count*")
count_dat <- map_df(count_f, function(x) 
  read_csv(sprintf("%s/%s", count_path, x), col_types="cd")) %>% 
  bind_rows() 
count_dat2 <- count_dat %>% filter(!is.na(sample_acc))

stopifnot(nrow(count_dat2)==length(unique(count_dat2$sample_acc)))

nrow(count_dat2)
count_dat2 %>% filter(read_count > 100000, !is.na(read_count)) %>% nrow()

count_dat2 %>% write_csv(sprintf("data/02_metadata/qc/all_%s_read_counts.csv", organism))

