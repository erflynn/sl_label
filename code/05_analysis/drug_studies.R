require('tidyverse')

drugbank_dat <- read_tsv("../labeling/geo2drug/data/02_labeled_data/drugbank_mapped_gse.txt")
no_cl_dat <- read_tsv("../labeling/geo2drug/data/02_labeled_data/non_cell_line_gse.txt")

load("data/summary_files.RData")
drugbank_dat %>% head()
drugbank_dat <- drugbank_dat %>% mutate(class=substr(ATC, 1,1))


# human hormone
hormone <- drugbank_dat %>% filter(class=="G" & organism=="human") %>% select(dbID, name, ATC, gse, gpl) %>% 
  inner_join(h_summ %>% filter(labeling_method=="exprsex") %>% select(study, sex) %>% rename(study_sex=sex), by=c("gse"="study"))

cancer_studies <- drugbank_dat %>% filter(class=="L" & organism=="human") %>% select(dbID, name, ATC, gse, gpl) %>% 
  inner_join(h_summ %>% filter(labeling_method=="exprsex") %>% select(study, sex) %>% rename(study_sex=sex), by=c("gse"="study"))
cancer_no_cl <- cancer_studies %>% left_join(no_cl_dat %>% rename(gse=gse_to_keep))
table(cancer_no_cl$study_sex)

neuro_hum <- drugbank_dat %>% filter(class=="N" & organism=="human") %>% select(dbID, name, ATC, gse, gpl) %>% 
  inner_join(h_summ %>% filter(labeling_method=="exprsex") %>% select(study, sex, num_samples) %>% rename(study_sex=sex), by=c("gse"="study"))
neuro_hum %>% write_csv("data/human_neur.csv")
neuro_hum %>% left_join(no_cl_dat %>% rename(gse=gse_to_keep))

# mouse neuro
neuro <- drugbank_dat %>% filter(class=="N" & organism=="mouse")
neur_w_summ <- neuro %>% select(dbID, name, ATC, gse, gpl) %>% 
  inner_join(m_summ %>% filter(labeling_method=="exprsex") %>% select(study, sex, num_samples) %>% rename(study_sex=sex), by=c("gse"="study"))
neur_w_summ %>% write_csv("data/mouse_neur.csv")

cardiac <- drugbank_dat %>% filter(class=="C" & organism=="mouse") %>% select(dbID, name, ATC, gse, gpl) %>% 
  inner_join(m_summ %>% filter(labeling_method=="exprsex") %>% select(study, sex) %>% rename(study_sex=sex), by=c("gse"="study"))
cardiac %>% write_csv("data/mouse_cardiac.csv")
