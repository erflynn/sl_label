# 02a_drug_study_labels.R
# E Flynn
# Updated 8/3/2020
# Label studies with drug mentions.
# 
# TODO:
# - it appears that these counts are LOWER than expected
# - I think this may be due to how we are grabbing metadata

require('tidyverse')
require('tidytext')
source("code/01_metadata/03_map/00_mapping_utils.R")
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
  mutate(across(c(title,description), ~ifelse(is.na(.) | . == "No description.", "", .))) %>%
  mutate(str=paste(title, description, sep=" ")) 

comb_data2 <- comb_data %>% 
  select(organism, data_type, study_acc, str) %>%
  unique()

# ----- drug data ------ #
# create the synonym DF
drug_full_info <- read.delim("../labeling/geo2drug/data/00_db_data/drugbank_parsed.txt")
drug_info_df <- read.delim("../labeling/geo2drug/data/00_db_data/drugbank_vocab_no_nutra.txt")

# only match drug data with more than 3 characters and 1 word 
#  TODO - what about drugs that are phrases?
drug_info_df$numchar <- sapply(drug_info_df$name, nchar)
short_names <- filter(drug_info_df, numchar <= 3)
names_to_match <- filter(drug_info_df, numchar > 3)
drug_info_df$nwords <- sapply(drug_info_df$name, function(x) length(strsplit(x, " ")[[1]]))

study_data <- comb_data2 %>% mutate(src_col="study")
study_lab <- labelNgram(study_data, drug_info_df)

#  problem - discarding aa from the drugbank classification
drug_stopwords <- c("same", "dmso", "water", "sage", "camp", "balance",  
                                      "biotin", "ethanol", "cyclo",
                                      "beam", "carbon", "silica", "amino acids", "ozone",
                                      "glucose", "d-glucose", "oxygen", "nitrogen", "l-glutamine", "sucrose")

jake_stopwords <- read.delim("../labeling/geo2drug/data/00_db_data/jake_stopwords.txt", head=FALSE)$V1


study_lab2 <- study_lab %>% filter( ! name  %in% union(drug_stopwords, jake_stopwords))
study_mapped <- study_lab2 %>% left_join(study_data %>% select(-src_col), by="str")

length(unique(study_mapped$study_acc)) # 7665
length(unique(study_mapped$dbID)) # 1104 


# join with ATC
drug_data_study <- study_mapped %>%
  select(organism, data_type, study_acc, dbID) %>%
  unique() %>%
  inner_join(drug_full_info %>% select(dbID, name, ATC), by="dbID") 

length(unique(drug_data_study$dbID))
length(unique(drug_data_study$study_acc)) 
drug_data_study %>% group_by(organism, data_type) %>% count()

# note - this is long format (each row is a "study, drug" not a "study)
drug_data_study %>% 
  arrange(organism, data_type, study_acc, name) %>%
  write_csv("data/02_labeled_data/study_to_drug.csv")



# --- counts + assessment for paper --- #
length(unique(comb_data2$study_acc)) # 44184
head(drug_data_study)
length(unique(drug_data_study$study_acc)) # 7665

non_drug <- comb_data %>%
  select(-str) %>%
  anti_join(drug_data_study, by="study_acc") %>%
  group_by(study_acc) %>%
  mutate(organism=paste(unique(organism), collapse=";"),
            data_type=paste(unique(data_type), collaspe=";"),
            title=paste(unique(title), collapse=";"),
            description=paste(unique(description), collapse=";")) %>%
  unique() %>%
  ungroup()

drug <- comb_data %>%
  inner_join(drug_data_study %>% select(-ATC)) %>%
  group_by(study_acc) %>%
  mutate(organism=paste(unique(organism), collapse=";"),
        data_type=paste(unique(data_type), collaspe=";"),
        dbID=paste(unique(dbID), collapse=";"),
        name=paste(unique(name), collapse=";"),
        title=paste(unique(title), collapse=";"),
        description=paste(unique(description), collapse=";")) %>%
  unique() %>% 
  ungroup()

comb_data %>%
  inner_join(drug_data_study %>% select(-ATC)) %>%
  distinct(dbID) %>% nrow()

set.seed(727)
non_drug %>% ungroup() %>% sample_n(200) %>% write_csv("data/eval_data/non_drug_study.csv")
drug %>% ungroup() %>% sample_n(200) %>% write_csv("data/eval_data/drug_study.csv")