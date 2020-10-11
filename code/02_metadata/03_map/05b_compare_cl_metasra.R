# // TODO: fix so not using comb_metadata to map samples to data types, etc

require('rjson')
require('tidyverse')

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", 
                          col_types="cccccccdld") 

mapped_to_cl <- read_csv("data/02_labeled_data/cell_line_mapping.csv")
cvcl_map <- fromJSON(file="../metasra_cvcl_mappings.json") 

cvcl_df <- do.call(rbind, lapply(cvcl_map, function(x) 
  paste(x$mapped_terms, collapse=";")))

load("data/10_comparison/metasra_data.RData") # --> mapped_df

# // TODO: what is the difference btw cl_mapped2 and cvcl_df???
cell_mapped <- mapped_df %>% 
  filter(str_detect(term_id, "CVCL")) %>%
  mutate(term_id2=str_replace_all(tolower(term_id), ":", "_"))
cl_mapped2 <- cell_mapped %>% 
  left_join(samp_to_runs, by=c("sample_accession"="sample")) %>%
  dplyr::rename(acc=run) 
head(cl_mapped2)

human_rnaseq <- comb_metadata %>% filter(data_type=="rnaseq", organism=="human", num_reads > 100000, present)
mapped_human_rnaseq <- mapped_to_cl %>% semi_join(human_rnaseq)
cl_mapped3 <- cl_mapped2 %>% semi_join(human_rnaseq, by=c("acc"="sample_acc"))
cl_mapped3 %>% anti_join(mapped_human_rnaseq, by=c("acc"="sample_acc")) %>%
  sample_n(10)
length(setdiff(cl_mapped3$acc, mapped_human_rnaseq$sample_acc)) # 6611
length(setdiff(mapped_human_rnaseq$sample_acc, cl_mapped3$acc)) # 15459

both_cl <- cl_mapped2 %>% inner_join(mapped_human_rnaseq, by=c("acc"="sample_acc"))
both_cl %>% filter(term_id2==accession) %>% nrow() # 14390/18821 (76.5%)
multi_map <- both_cl %>% filter(term_id2!=accession & str_detect(accession, ";")) 
some_overlap <- multi_map %>% filter(str_detect(accession,term_id2)) # 629
no_overlap <- multi_map %>% filter(!str_detect(accession,term_id2)) # 326


# metrics:
#  [x] compare to previous
#  [x] fraction of cell line fields mapping
#  [x] compare to metaSRA
#  - 100 randomly selected that do not map but contain a cell word 
#  - 100 randomly selected that do map