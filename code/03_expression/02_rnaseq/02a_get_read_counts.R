# read in the .csv file and grab the read counts by summing the columns

require('tidyverse')
args <- commandArgs(trailingOnly=TRUE)
organism <- args[1]
my_idx <- as.numeric(args[2])

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}

my_f <- str_replace_all(list.files(path=sprintf("data/rnaseq/%s/01_count_mat", organism), 
                                   pattern="^[E|D|S]RP*"), ".csv", "")
studies <- extractChunk(my_f, my_idx, SIZE.CHUNK=500)

get_counts <- function(study_id){
  dat <- read_csv(sprintf("%s.csv", study_id))
  summed_cols <- colSums(dat %>% select(-GENEID), na.rm=T)
  my_df <- tibble("sample_acc"=names(summed_cols),  
                  "study_acc"=rep(study_id, length(summed_cols)), 
                  "count"=summed_cols)
  return(my_df)
}

count_df <- studies %>% map_df(get_counts) %>% bind_rows()
count_df %>% 
  write_csv(sprintf("data/rnaseq/%s/02_read_counts/count_%s.csv", organism, my_idx))
