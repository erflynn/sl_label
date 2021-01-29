# drug_gse_labeling.R
# E Flynn
# 4/14/2019

require('tidyverse')
require('tidytext')
options(stringsAsFactors = FALSE)

# read in the metadata and clean up missing fields
gse_data <- read.csv("data/01_metadata/human_experiment_metadata.csv", stringsAsFactors = FALSE) %>% 
  as_tibble() %>% 
  select(-date) %>%
  mutate(organism="human") %>%
  bind_rows(read.csv("data/01_metadata/mouse_experiment_metadata.csv", stringsAsFactors = FALSE) %>% 
              as_tibble() %>% 
              select(-date) %>% 
              mutate(organism="mouse"))
  
sra_data <- read_csv("data/01_metadata/human_rnaseq_experiment_metadata.csv") %>% 
  select(-date) %>%
  mutate(organism="human") %>%
  bind_rows(read_csv("data/01_metadata/mouse_rnaseq_experiment_metadata.csv") %>% 
              select(-date) %>%
              mutate(organism="mouse"))

# remove SRA data from compendia and create a string field
comb_data <- gse_data %>% anti_join(sra_data, by="study_acc") %>%
  bind_rows(sra_data) %>%
  mutate(across(c(title,description), ~ifelse(is.na(.) | . == "No description.", "", .))) %>%
  mutate(str=paste(title, description, sep=" ")) 

comb_data2 <- comb_data %>% select(organism, study_acc, str)

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

# fix tabs and trailing punctuation
text_df <- gse_data[,c("gse", "str")]
text_df$str <- sapply(text_df$str, function(x)
  {y <- gsub(";\t", " ", x);  ## fix ;\t
  z <- gsub(' [[:punct:]]|[[:punct:]] ', ' ', y);  ##remove trailing punctuation
  return(z)})

# separate out the gse data into unigrams
gse_unigrams <- text_df %>% unnest_tokens(word, str, token=stringr::str_split, pattern = " ") 
gse_unigrams2 <- filter(gse_unigrams, str_length(word) > 3 & is.na(as.numeric(word))) # remove short words or numbers 

# map with a join
comb_names <- inner_join(filter(drug_info_df, nwords==1), gse_unigrams2, by=c("name"="word")) %>% distinct()
length(unique(comb_names$gse)) 
drugs <- unique(comb_names$name ) 
gse_map_counts <- gse_unigrams %>% filter (word %in% drugs) %>% group_by(word) %>% summarise(total=n())  %>% arrange(desc(total))


# try some bigrams
bigrams <- gse_data[,c("gse", "str")] %>% unnest_tokens(bigram, str, token = "ngrams",  n = 2)
trigrams <- gse_data[,c("gse", "str")] %>% unnest_tokens(trigram, str, token = "ngrams", n = 3)

comb_names_bi <- inner_join(filter(drug_info_df, nwords==2), bigrams, by=c("name"="bigram")) %>% distinct()
comb_names_tri <- inner_join(filter(drug_info_df, nwords==3), trigrams, by=c("name"="trigram")) %>% distinct()

#  problem - discarding aa from the drugbank classification
drug_stopwords <- c("same", "dmso", "water", "sage", "camp", "balance",  "biotin")
# need to keep in ATRA, SAHA

# pull in additional curated stopwords
jake_stopwords <- read.delim("data/00_db_data/jake_stopwords.txt", head=FALSE)$V1
comb_names2 <- filter(select(comb_names, c("name", "dbID", "gse")), ! name  %in% union(drug_stopwords, jake_stopwords))
length(unique(comb_names2$gse)) # 9735
length(unique(comb_names2$dbID)) # 1104


# join with ATC
comb_names3 <- left_join(comb_names2, drug_full_info[,c("dbID", "ATC")], by="dbID")

# NOW TO DISAMBIGUATE - currently has synonymns as names
drug_data_gse <- inner_join(drug_full_info[,c("dbID", "name")], select(comb_names3, c("dbID", "ATC", "gse")))


drug2 <- filter(drug_data_gse, ! name %in% 
                  c("Glucose", "D-glucose", "Oxygen", "Nitrogen", "L-Glutamine", "Sucrose"))
length(unique(drug2$dbID)) # [1] 1098
length(unique(drug2$gse)) # [1] 8834

# add a mouse or human flag
drug2_gse <- left_join(drug2, gse_data %>% select(gse, gpl, organism, study_type))

# FYI these are multimapping!
write.table(drug2_gse, file="data/02_labeled_data/drugbank_mapped_gse.txt", row.names=FALSE, sep="\t")

d2 <- drug2_gse %>% select(gse, organism, study_type) %>% unique()
d2 %>% group_by(organism, study_type) %>% count()
# human    oligo       3883
# human    seq         1017
# mouse    oligo       2862
# mouse    seq         1072

# write out download lists
drug2_gse %>% filter(organism=="mouse" & study_type=="oligo") %>% select(gse) %>% unique() %>%
  write_csv("data/01_sample_lists/mouse_oligo_drug.csv") # 2862
drug2_gse %>% filter(organism=="human" & study_type=="oligo") %>% select(gse) %>% unique() %>% 
  write_csv("data/01_sample_lists/human_oligo_drug.csv") # 3883


