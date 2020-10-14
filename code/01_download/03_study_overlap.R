# goal:
# - any two samples that share a study are listed as the same study group
# results in a file "data/01_sample_lists/sample_to_study_grp.csv" that contains this info

require('tidyverse')

mapping <- read_csv("data/01_sample_lists/sample_to_study.csv") 
# distinct row: sample, study, organism

# count by sample
counts_per_sample <- mapping %>% 
  filter(sample_acc != "") %>%
  group_by(organism, sample_acc) %>% count()

# get a list of the samples that map to more than one study
multi_sample <- counts_per_sample %>% 
  filter(n > 1) %>% left_join(mapping) 

# get a list of the studies that overlap for these
mult_samp_to_study <- multi_sample %>% 
  select(-n, -data_type, -src) %>%
  group_by(organism, sample_acc) %>% 
  arrange(study_acc) %>%
  summarize(study_acc=paste(unique(study_acc), collapse=";")) %>%
  ungroup() %>%
  select(organism, study_acc) %>% 
  unique() %>% 
  arrange(organism, study_acc) %>%
  mutate(ord=1:n())

agg_studies <- mult_samp_to_study %>% 
  separate_rows(study_acc) %>% 
  group_by(organism, study_acc) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  filter(n>5) # 82 studies are in more than 5 groups



study_to_id <- mult_samp_to_study %>% 
  separate_rows(study_acc) %>% 
  anti_join(agg_studies) %>%
  group_by(organism, study_acc) %>% 
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
new_df <- data.frame(cbind("study"=unlist(names(study_to_ord)), 
                           "new_id"=unlist(study_to_ord)))
ord_map <- data.frame(cbind("new_id"=unlist(names(new_ord)), 
                            "final_id"=unlist(new_ord)))
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

still_multi <- new_mapping2 %>% group_by(organism, sample_acc) %>% filter(n() >1) 
stopifnot(nrow(still_multi)==0)



# one line is an organism, sample
stopifnot(nrow(new_mapping2)==(new_mapping2 %>% distinct(organism, sample_acc) %>% nrow()))


# the result of this is that we have mapped each sample to a study_str
#  which is a group of studies with a distinct set of samples
new_mapping2 %>% rename("study_grp"=study_str) %>% 
  write_csv("data/01_sample_lists/sample_to_study_grp.csv")