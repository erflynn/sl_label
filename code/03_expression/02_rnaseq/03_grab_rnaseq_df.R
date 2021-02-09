require('tidyverse')
require('glmnet')
require('data.table')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
sample_f <- args[2]
outfile <- args[3]

grab_samples <- function(study_id){
  sample_list2 <- sample_list %>% filter(study_acc==study_id)
  my.f <- sprintf("data/03_expression/rnaseq/%s/01_study_mat/%s.csv", 
                  prefix, study_id)
  if (!file.exists(my.f) | file.info(my.f)$size==0){
    return(NA)
  }
  # what to do if the file doesn't exist
  dat <- fread(my.f, data.table=FALSE)
  overlap.s <- intersect(colnames(dat), sample_list2$sample_acc)
  return(dat[,c("gene_name", overlap.s)])
}

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}

#all_dfs <- lapply(head(unique(sl_sm$study_acc), 10), grab_samples)
#all_dfs2 <- all_dfs[!is.na(all_dfs)]

# read in a file with the list of samples to load
sample_list <- fread(sample_f, data.table=FALSE)


list_studies <- unique(sample_list$study_acc)
num_chunks <- ceiling(length(list_studies)/20)
all_df <- data.frame()

# iterate through and full join
for (i in 1:num_chunks){
  print(i)
  sm_list_studies <- extractChunk(list_studies, i-1, 20)
  dfs <- lapply(sm_list_studies, grab_samples)
  dfs2 <- dfs[!is.na(dfs)]
  chunk_df <- dfs2 %>% reduce(full_join, by="gene_name")
  if (!is.na(chunk_df)){
    save(chunk_df, file=sprintf("data/03_expression/rnaseq/%s/04_df/%s_%s.RData", prefix, outfile, i))
    if (ncol(all_df)==0){
      all_df <- data.frame(chunk_df)
    } else {
      all_df <- full_join(all_df, data.frame(chunk_df), by="gene_name")
    }
  }
}

save(all_df, 
     file=sprintf("data/03_expression/rnaseq/%s/04_df/%s_rnaseq_%s.RData", prefix, prefix, outfile))






