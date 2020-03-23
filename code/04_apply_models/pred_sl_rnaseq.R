require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
idx <- as.numeric(args[2])

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}


slStudy <- function(study_id){
  df <- read_csv(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id))  %>% as.data.frame()
  
  # rotate
  
  # select genes
  
  # label
  
  return() # return labels
}

  
  