require('tidyverse')
require('glmnet')
require('data.table')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

sl <- read_csv(sprintf("data/%s/02_sample_lists/%s_rnaseq_sex_lab.csv", 
                       prefix, prefix))
mapping2 <- read_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample2.csv", prefix))
mapping3 <- mapping2 %>% filter(present)
sl2 <- sl %>% filter(sample_acc %in% mapping3$sample_acc)
sl3 <- sl2 %>% group_by(study_acc) %>% filter(n() >=5) %>% sample_n(5)
sl3.2 <- sl2 %>% group_by(study_acc) %>% filter(n() <5) 
sl4 <- rbind(sl3, sl3.2)

sl_sm <- sl4

grab_samples <- function(study_id){
  sample_list <- sl_sm %>% filter(study_acc==study_id)
  my.f <- sprintf("data/rnaseq/%s/01_study_mat/%s.csv", 
                  prefix, study_id)
  if (!file.exists(my.f) | file.info(my.f)$size==0){
    return(NA)
  }
  # what to do if the file doesn't exist
  dat <- fread(my.f, data.table=FALSE)
  overlap.s <- intersect(colnames(dat), sample_list$sample_acc)
  return(dat[,c("gene_name", overlap.s)])
}

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}

#all_dfs <- lapply(head(unique(sl_sm$study_acc), 10), grab_samples)
#all_dfs2 <- all_dfs[!is.na(all_dfs)]
list_studies <- unique(sl_sm$study_acc)
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
    if (ncol(all_df)==0){
      all_df <- data.frame(chunk_df)
    } else {
      all_df <- full_join(all_df, data.frame(chunk_df), by="gene_name")
    }
  }
}

save(all_df, 
     file=sprintf("data/%s/04_sl_input/%s_rnaseq_sl_in2.RData", prefix, prefix))






