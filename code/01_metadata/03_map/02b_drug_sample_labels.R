# 02b_drug_sample_labels.R
# E Flynn
# Updated 8/3/2020
#
# Label samples that have drug info based on the refine-bio treatment field

require('tidyverse')
require('tidytext')
source("code/01_metadata/03_map/00_mapping_utils.R")

options(stringsAsFactors = FALSE)

# read in the sample metadata
compendia_metadata <- read.csv("data/01_metadata/human_metadata.csv") %>%
  mutate(organism="human", data_type="microarray") %>%
  bind_rows(read.csv("data/01_metadata/mouse_metadata.csv") %>%
              mutate(organism="mouse", data_type="microarray") ) %>%
  select(organism, data_type, acc, compound, trt, title)
rnaseq_metadata <- read.csv("data/01_metadata/human_rnaseq_sample_metadata.csv") %>%
  mutate(organism="human", data_type="rnaseq") %>%
  bind_rows(read.csv("data/01_metadata/mouse_rnaseq_sample_metadata.csv") %>%
              mutate(organism="mouse", data_type="rnaseq") ) %>%
  select(organism, data_type, acc, compound, trt, title)

# put together, removing the rnaseq samples from compendia
comb_trt_data <- compendia_metadata %>%
  anti_join(rnaseq_metadata, by="acc") %>%
  bind_rows(rnaseq_metadata) %>%
  rename(sample_acc=acc) %>%
  as_tibble() %>%
  mutate(across(c(compound, trt, title), tolower))

trt_words <- comb_trt_data %>% group_by(trt) %>% dplyr::count()
trt_words2 <- trt_words %>% filter(trt != "")

# -- missing a lot of IL
trt <- comb_trt_data %>% filter(trt != "" & 
                      is.na(as.numeric(trt)) & 
                      nchar(trt)>3) 
trt_unique <- trt %>% group_by(trt) %>% count() %>% arrange(desc(n))

# //TODO: commenting out b/c unclear, come back to!
#mapping <- read.csv("data/01_metadata/human_exp_to_sample.csv")
#drug_mapped <- read_csv("data/human_sample_drug_mapped.csv") #// TODO: where is this?!!!
#trt2 <- trt %>% left_join(mapping, by=c("acc"="sample_acc"))
#trt3 <- trt2 %>% left_join(drug_mapped, by=c("acc"="gsm"))
#trt3 %>% filter(!is.na(dbID)) %>% select(study_acc)  %>% unique() %>% nrow()

# ----- drug data ------ #
# create the synonym DF
drug_full_info <- read.delim("../labeling/geo2drug/data/00_db_data/drugbank_parsed.txt")
drug_info_df <- read.delim("../labeling/geo2drug/data/00_db_data/drugbank_vocab_no_nutra.txt")

# read in the
ctl_vocab <- c("none", "no", "control", "untreated", "dmso", "na", "placebo", "saline", "pbs", "mock",
               "baseline", "unstimulated", "etoh", "ethanol", "ctrl", "non-treated", "vehicle", "ctl")

ctl_dat <- cbind("name"=ctl_vocab, "dbID"=rep("ctl", length(ctl_vocab) )) %>% as_tibble()
non_drug_vocab <- c("reprogramming", "stimulation", "injured", "radiation", "transdifferentiation", 
                    "injury", "heatshock", "heat")
non_drug <- cbind("name"=non_drug_vocab, "dbID"=rep("non-drug", length(non_drug_vocab))) %>% as_tibble()

ctl_dat2 <- ctl_dat %>% bind_rows(non_drug) %>% mutate(nwords=1) 
  

# only match drug data with more than 3 characters and 1 word 
#  TODO - what about drugs that are phrases?
drug_info_df$numchar <- sapply(drug_info_df$name, nchar)
short_names <- filter(drug_info_df, numchar <= 3)
names_to_match <- filter(drug_info_df, numchar > 3)
drug_info_df$nwords <- sapply(drug_info_df$name, function(x) length(strsplit(x, " ")[[1]]))

trt2 <- trt %>%
  pivot_longer(compound:title, names_to="src_col", values_to="str")

ctl_names <- labelNgram(trt2, ctl_dat2, remove_short = FALSE)
comb_names <- labelNgram(trt2, drug_info_df)

comb_names %>% select(dbID, str) %>% unique() %>%
  group_by(str) %>% count() %>% filter(n>1) %>% nrow() # 390 multimap

drug_stopwords <- c("same", "dmso", "water", "sage", "camp", "balance",  "biotin", "ethanol", "cyclo",
                    "beam", "carbon", "silica", "amino acids", "ozone",
                    "glucose", "d-glucose", "oxygen", "nitrogen", "l-glutamine", "sucrose") # need to keep in ATRA, SAHA
jake_stopwords <- read.delim("../labeling/geo2drug/data/00_db_data/jake_stopwords.txt", head=FALSE)$V1

mapped_df <- trt2 %>% 
  inner_join(comb_names, by=c("src_col", "str")) %>%
  filter(! name  %in% union(drug_stopwords, jake_stopwords)) 
length(unique(mapped_df$dbID)) # 436 drugs 

# TODO: map some 3char drugs? e.g. il1, il7, ...

ctl_mapped <- trt2 %>%
  inner_join(ctl_names, by=c("src_col", "str")) 

ctl_mapped %>% filter(dbID=="ctl", src_col=="trt") %>% distinct(sample_acc) %>% nrow() # 31155
ctl_mapped %>% filter(dbID=="non-drug" & src_col=="trt") %>% distinct(sample_acc) %>% nrow() # 7262

# how many have control -AND- drug
mapped_df %>% inner_join(ctl_mapped %>% 
                           filter(dbID=="ctl") %>%
                           select(sample_acc, src_col, name), 
                         by=c("sample_acc", "src_col")) %>% 
  #filter(src_col=="trt") %>%
  select(src_col, str, name.x, name.y)
# 1068 total, 600 in treatment

#mapped_df %>%
#  group_by(name) %>% 
#  count() %>% 
#  arrange(desc(n)) %>%
#  View()

# add ATC label and convert "name" (which could be a synonym) to the true name
drug_data_sample <- mapped_df %>%
  select(organism, data_type, sample_acc, dbID) %>%
  unique() %>%
  inner_join(drug_full_info %>% select(dbID, name, ATC), by="dbID") 

length(unique(drug_data_sample$dbID)) # 442
length(unique(drug_data_sample$sample_acc)) # --> 13653

# how many with multiple drugs?
drug_data_sample %>% 
  group_by(sample_acc) %>% 
  count() %>% 
  arrange(desc(n)) %>%
  filter(n > 1) %>%
  nrow() # 783 mapped to multiple drugs!

ggplot(drug_data_sample %>% 
         group_by(sample_acc) %>% 
         count(), 
       aes(x=n)) + 
  geom_histogram(stat="count")

drug_data_sample %>% write_csv(sprintf("data/rb_sample_drug_mapped.csv"))

mapped_all <- mapped_df %>%
  select(organism, data_type, sample_acc,  src_col, str, name, dbID ) %>%
  mutate(type="drug") %>%
  bind_rows(ctl_mapped %>% select(organism, data_type, sample_acc,  src_col, str, name, dbID ) %>%
              mutate(type=dbID)) 
mapped_all %>% write_csv("data/all_drug_ctl_annot_sample.csv")

mapped_lab <- mapped_all %>% 
  group_by(organism, data_type, sample_acc) %>%
  mutate(sample_type=case_when(
    any(type=="non-drug") ~ "non-drug",
    any(type=="ctl") & any(type=="drug") ~ "both",
    all(type=="drug") ~ "drug",
    all(type=="ctl") ~ "ctl"
  )) %>% 
  select(-src_col, -str, -name, -dbID, -type) %>%
  unique() 

table(mapped_lab$sample_type)

mapped_lab %>% write_csv("data/sample_mapped_drug_type.csv")

# --------- end of sample labeling ----------- #


# ---- what is the overlap with study data ? ---- #
sample_data <- read_csv("data/rb_sample_drug_mapped.csv")
study_data <- read_csv("data/02_labeled_data/study_to_drug.csv")

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")

sample_study <- comb_metadata %>% 
  select(sample_acc, study_acc) %>% 
  separate_rows(study_acc,  sep=";")

sample_data2 <- sample_data %>% left_join(sample_study, by="sample_acc")
sample_studies <- unique(sample_data2$study_acc) # 1095
study_mentions <- unique(study_data$study_acc) # 7665
length(intersect(sample_studies, study_mentions)) # 757
length(setdiff(sample_studies, study_mentions)) # 338
length(setdiff(study_mentions, sample_studies)) # 6908

# --- are they the same drug? --- #
samp_study_overlap <- sample_data2 %>% 
  inner_join(study_data %>% select(-ATC) , by=c("organism", "data_type", "study_acc"))

table(samp_study_overlap$dbID.x==samp_study_overlap$dbID.y) # 5644 false, 9927 true

samp_study_compare <- samp_study_overlap %>%
  group_by(organism, data_type, study_acc) %>%
  mutate(dbIDx=paste(unique(dbID.x), collapse=";"),
         dbIDy=paste(unique(dbID.y),collapse=";"),
         dbID_overlap=length(intersect(dbID.x, dbID.y)),
         dbID_union=length(union(dbID.x, dbID.y))) %>%
  select(organism, data_type, study_acc, dbIDx, dbIDy, dbID_overlap, dbID_union) %>%
  unique()

samp_study_compare %>% filter(dbID_overlap==dbID_union) # 557 studies have the same drugbank annotations

samp_study_compare %>% filter(dbID_overlap==0) %>% nrow() # 42 have no overlap in the drug annotations
sample_data2 %>%
  semi_join(samp_study_compare %>% filter(dbID_overlap==0)) %>%
  select(-sample_acc) %>% 
  unique() %>%
  inner_join(study_data %>% select(-organism, -data_type), by=c("study_acc")) %>%
  View()

samp_study_compare %>%
  filter(dbID_overlap!=dbID_union & dbID_overlap > 0) %>%
  arrange(desc(dbID_union), dbID_overlap) # 159 overlap by at least one drug but not all


# ---------- set up assessment data ------ #
nrow(trt) # number of samples with treatment fields
trt_unique # number of unique

# just control terms
mapped_all %>% 
  filter(type=="ctl" & src_col=="trt") %>% 
  anti_join(mapped_all %>% filter(type=="drug" & src_col=="trt"), by="sample_acc") %>%
  distinct(sample_acc) %>% nrow()
mapped_all %>% 
  filter(type=="ctl" & src_col=="trt") %>% 
  anti_join(mapped_all %>% filter(type=="drug" & src_col=="trt"), by="sample_acc") %>%
  distinct(str) %>% nrow()

# drugs and control terms
mapped_all %>% 
  filter(type=="ctl" & src_col=="trt") %>% 
  semi_join(mapped_all %>% filter(type=="drug" & src_col=="trt"), by="sample_acc") %>%
  distinct(sample_acc) %>% nrow()
mapped_all %>% 
  filter(type=="ctl" & src_col=="trt") %>% 
  semi_join(mapped_all %>% filter(type=="drug" & src_col=="trt"), by="sample_acc") %>%
  distinct(str) %>% nrow()

# just drugs
mapped_all %>% 
  filter(type=="drug" & src_col=="trt") %>% 
  anti_join(mapped_all %>% filter(type=="ctl" & src_col=="trt"), by="sample_acc") %>%
  distinct(sample_acc) %>% nrow()

mapped_all %>% 
  filter(type=="drug" & src_col=="trt") %>% 
  anti_join(mapped_all %>% filter(type=="ctl" & src_col=="trt"), by="sample_acc") %>%
  distinct(str) %>% nrow()

mapped_all %>% 
  filter(type=="drug" & src_col=="trt") %>% 
  anti_join(mapped_all %>% filter(type=="ctl" & src_col=="trt"), by="sample_acc") %>%
  distinct(dbID) %>% nrow()

unmapped <- trt %>%
  anti_join(
    mapped_all %>% 
      filter(type %in% c("drug", "ctl") & 
               src_col=="trt"), 
    by="sample_acc"
  )

unmapped %>% distinct(trt) %>% nrow()
head(unmapped)

# additional mapping
unmapped %>% semi_join(
  mapped_all %>% 
    filter(type %in% c("drug") & 
             src_col %in% c("title", "compound")), 
  by="sample_acc"
)

# ----- assessment ----- #

ctl_samples <- trt %>% semi_join(
  mapped_all %>% filter(type=="ctl") ,
  by="sample_acc"
)
drug_samples <- trt %>% inner_join(
  mapped_all %>% 
    select(-src_col, -organism, -data_type, -str) %>% 
    unique() %>% 
    filter(type=="drug") ,
  by="sample_acc")
unmapped_samples <- trt %>% 
  anti_join(ctl_samples, by="sample_acc") %>% 
  anti_join(drug_samples, by="sample_acc")

set.seed(252)
ctl_samples %>% sample_n(200) %>% write_csv("data/eval_data/control_samples.csv")
drug_samples %>% sample_n(200) %>% write_csv("data/eval_data/drug_samples.csv")
unmapped_samples %>% sample_n(200) %>% write_csv("data/eval_data/unmapped_samples.csv")

#//TODO: maybe need to add study info to these to get a better assessment?
