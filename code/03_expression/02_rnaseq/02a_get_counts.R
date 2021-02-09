# read in all the data and get counts

library(tidyverse)
library(vroom)
args <- commandArgs(trailingOnly=TRUE)

organism <- args[1]
idx <- as.numeric(args[2])

setwd(sprintf("data/03_expression/rnaseq/%s/01_study_mat/", organism))
my_f <- list.files(pattern="*.csv")

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}
grab_chunk <- function(idx){
  my_f0 <- extractChunk(my_f, idx, SIZE.CHUNK = 500)
  my_f_r <- lapply(my_f0, function(f) {
    df <- vroom(f);
    if (ncol(df)==0){
      return(NA)
    }
    df %>% filter(gene_name=="read_count") %>% select(-gene_name)})
  
  df2 <- do.call(cbind,my_f_r)
  counts <- unlist(df2[1,])
  cols <- colnames(df2)
  tb <- tibble("sample_acc"=cols, "read_count"=counts)
  return(tb)
}

my_chunk <- grab_chunk(idx)
my_chunk %>% write_csv(sprintf("../02_counts/count_chunk_%s.csv", idx))