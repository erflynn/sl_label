# put together the study metadata

require('tidyverse')
require('tidytext')
#source("code/01_metadata/03_map/00_mapping_utils.R")
options(stringsAsFactors = FALSE)

# read in the metadata and clean up missing fields
compendia_data <- read.csv("data/01_metadata/human_experiment_metadata.csv") %>% 
  as_tibble() %>% 
  select(-date) %>%
  mutate(organism="human") %>%
  bind_rows(read.csv("data/01_metadata/mouse_experiment_metadata.csv") %>% 
              as_tibble() %>% 
              select(-date) %>% 
              mutate(organism="mouse")) %>%
  mutate(data_type="microarray") %>%
  unique()

sra_data <- read_csv("data/01_metadata/human_rnaseq_experiment_metadata.csv") %>% 
  select(-date) %>%
  mutate(organism="human") %>%
  bind_rows(read_csv("data/01_metadata/mouse_rnaseq_experiment_metadata.csv") %>% 
              select(-date) %>%
              mutate(organism="mouse")) %>%
  mutate(data_type="rnaseq")

# remove SRA data from compendia and create a string field
comb_data <- compendia_data %>% 
  anti_join(sra_data, by="study_acc") %>%
  bind_rows(sra_data) %>%
  mutate(across(c(title,description), 
                ~ifelse(is.na(.) | . == "No description.", "", .))) %>%
  mutate(str=paste(title, description, sep=" ")) 

comb_data2 <- comb_data %>% 
  select(organism, data_type, study_acc, str) %>%
  unique()
# what study data is present

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", col_types="cccccccdld")
distinct_studies <- comb_metadata %>% 
  separate_rows(study_acc, sep=";") %>% 
  distinct(study_acc, data_type, organism)
# 33,470

# distinct present studies?
distinct_studies2 <- comb_metadata %>% 
  filter(data_type!="rnaseq" | (num_reads > 100000 & present))  %>% 
  separate_rows(study_acc, sep=";") %>% 
  distinct(study_acc, data_type, organism)
# 33,089

distinct_studies2 %>% group_by(organism, data_type) %>% count() 

comb_data3 <- comb_data %>% distinct(study_acc, title, description) # 44184

comb_data3 %>% inner_join(distinct_studies2) %>% nrow() # 33089
comb_data3 %>% anti_join(distinct_studies) # what happened to the 12k that aren't present...
comb_data3 %>% inner_join(distinct_studies2) %>% filter(title=="" | description=="") %>% nrow()
# not too many are missing data