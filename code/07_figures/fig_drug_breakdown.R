require("tidyverse")
drugbank_dat <- read_tsv("../labeling/geo2drug/data/02_labeled_data/drugbank_mapped_gse.txt")
no_cl_dat <- read_tsv("../labeling/geo2drug/data/02_labeled_data/non_cell_line_gse.txt")
human_dedup <- read_csv("data/01_metadata/human_compendia_dedup_mapping.csv")
refine_mapped <- read_csv("data/02_labeled_data/human_rnaseq_sample_cl.csv")
load("data/summary_files.RData")
drugbank_dat %>% head()

human_dedup2 <- human_dedup %>% separate_rows(study_str, sep=",")

drug_no_cl <- drugbank_dat %>%
  filter(organism=="human") %>%
  left_join(human_dedup2, by=c("gse"="study_str")) %>%
  anti_join(refine_mapped, by=c("sample_acc"="gsm"))

drug_no_cl %>% head()
drug_no_cl2 <- drug_no_cl %>% 
  select(gse, dbID, name, ATC) %>%
  unique() %>%
  left_join(h_summ %>% 
              filter(labeling_method=="exprsex"), by=c("gse"="study")) %>%
  select(gse, dbID, name, ATC, num_samples, num_f, num_m, sex) %>%
  filter(!is.na(sex))


drug_no_cl3 <- drug_no_cl2 %>% 
  mutate(class=substr(ATC, 1, 1)) %>%
  filter(!is.na(class)) %>%
  mutate(study_sex=case_when(
    num_f == 0 & num_m != 0 ~ "male-only",
    num_m == 0 & num_f != 0 ~ "female-only",
    TRUE ~ "mixed"
  )) %>%
  filter(num_samples > 10) 

counts_per_grp <-  drug_no_cl3 %>%
  
  group_by(class, study_sex) %>%
  count() 
  

ggplot(counts_per_grp, aes(x=factor(class), y=n, fill=class))+
  geom_histogram(stat="identity") +
  facet_wrap(. ~ study_sex)

study_metadata <- read.csv("data/01_metadata/human_experiment_metadata.csv")

drug_no_cl3 %>% filter(class %in% c( "L") & study_sex=="female-only") %>%
  left_join(study_metadata, by=c("gse"="study_acc")) %>%
  View()
