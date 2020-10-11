# study_sharing.R
#
# There are many samples in multiple studies. Because we want to prevent leakage 
# from training into testing, we need to make sure that all studies are completely 
# separate such that no two studies share a sample. This is a more stringent condition 
# than just simple deduplication and assignment to individual studies.
#
# To do this:
#  1) we remove "aggregate" studies which we define as studies that have overlap with >5 other studies
#  2) we group together all studies that share any samples and assign them a group study ID


require('tidyverse')
require('data.table')
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
data_type <- args[2]

# read in the mapping data
if (data_type == "compendia"){
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample.csv", prefix), data.table=FALSE)  
} else {
  MIN.READS <- 100000
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix), data.table=FALSE)  %>%
    filter(present & num_reads >= MIN.READS)
}



# count by sample
counts_per_sample <- mapping %>% 
  filter(sample_acc != "") %>%
  group_by(sample_acc) %>% count()

multi_sample <- counts_per_sample %>% 
  filter(n > 1) %>% left_join(mapping) %>% select(-n)

if (nrow(multi_sample)==0){
  print("no duplicates")
  stop()
}

mult_samp_to_study <- multi_sample %>% group_by(sample_acc) %>% 
  arrange(study_acc) %>%
  summarize(study_acc=paste(study_acc, collapse=",")) %>%
  ungroup() %>%
  select(study_acc) %>% 
  unique() %>% 
  arrange(study_acc) %>%
  mutate(ord=1:n())

agg_studies <- mult_samp_to_study %>% 
  separate_rows(study_acc) %>% 
  group_by(study_acc) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  filter(n>5) # 76 studies out of all of them, remove

study_to_id <- mult_samp_to_study %>% 
  separate_rows(study_acc) %>% 
  anti_join(agg_studies) %>%
  group_by(study_acc) %>% 
  mutate(ord=as.numeric(ord)) %>%
  arrange(ord) %>%
  summarize(ord_str=paste(ord, collapse=","),
            first_ord=unlist(ord)[1])
  
study_to_id2 <- study_to_id %>% 
  arrange(first_ord) %>% 
  select(-first_ord) %>% 
  rename(ord=ord_str)

# iterate through and assign each
study_to_ord = list()
new_ord=list()
for (i in 1:nrow(study_to_id2)){
  study_acc <- study_to_id2$study_acc[i]
  ord <- study_to_id2$ord[i]
  ords <- str_split(ord, pattern=",")[[1]]

  transformed_ords <- sapply(unlist(sapply(ords, function(x) {
    ifelse(x %in% names(new_ord), new_ord[x], x)})), as.numeric)
  
  names(transformed_ords) <- NULL
  min_ord=min(unlist(transformed_ords))
 
  for (j in 1:length(ords)){
    curr_ord=ords[j]
    new_ord[curr_ord]=min_ord
  }
  study_to_ord[study_acc] <- min_ord
}

# join together the two dfs
new_df <- data.frame(cbind("study"=unlist(names(study_to_ord)), "new_id"=unlist(study_to_ord)))
ord_map <- data.frame(cbind("new_id"=unlist(names(new_ord)), "final_id"=unlist(new_ord)))
new_df2 <- new_df %>% left_join(ord_map)
stopifnot(max((new_df2 %>% group_by(study) %>% count())$n)==1) # sanity check, the max should be one

# reformat 
new_df3 <- new_df2 %>%
  select(-new_id) %>%
  group_by(final_id) %>%
  arrange(final_id, study) %>%
  mutate(study_str=paste(study, collapse=",")) %>%
  rename(study_acc=study) %>%
  ungroup()

# add in the sample ID
new_mapping <- mapping %>%
  filter(sample_acc != "") %>%
  anti_join(agg_studies) %>%
  left_join(new_df3 %>% select(-final_id)) %>%
  mutate(study_str=ifelse(is.na(study_str), study_acc, study_str))

new_mapping2 <- new_mapping %>% select(-study_acc) %>% unique() 

still_multi <- new_mapping2 %>% group_by(sample_acc) %>% filter(n() >1) 
stopifnot(nrow(still_multi)==0)

# write out the new mapping
new_mapping2 %>% 
  write_csv(sprintf("data/01_metadata/%s_%s_dedup_mapping.csv", prefix, data_type))

