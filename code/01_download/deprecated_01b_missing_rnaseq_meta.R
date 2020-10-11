# 
# Code for grabbing lists of missing RNA-seq metadata
#
# TODO:
# - fix paths: they are relative to data/sra_out right now

require('tidyverse')
# ---- grab missing data ----- #
missing_f = list.files(pattern="missing*")
missing_dat <- do.call(rbind, lapply(missing_f, read_tsv, col_type="c", col_names=FALSE)) %>%
  as_tibble() # 9,966

m_dat2 <- missing_dat %>% 
  rename(ID="X1") %>%
  separate_rows(ID, sep="-") %>%
  arrange(ID) %>%
  unique() %>% #  9,601
  filter(str_detect(ID, "ER|SR|DR")) %>% # 9452
  mutate(type=case_when(
    str_detect(ID, "RX") ~ "experiment",
    str_detect(ID, "RR") ~ "run",
    str_detect(ID, "RS") ~ "sample"
  )) %>%
  filter(!is.na(type)) %>%
  mutate(id_len=nchar(ID))

# make sure the runs exist
all_runs <- read_csv("../rnaseq_runs_list/rnaseq_runs.csv", col_names=FALSE) %>%
  rename(ID=X1)
missing_runs <- m_dat2 %>% 
  filter(type=="run") %>% # 728
  semi_join(all_runs, by="ID") # --> 636
missing_runs %>% select(ID) %>% 
  write_csv("../missing_runs_to_download.csv", col_names=FALSE)

# samples
m_dat2 %>% 
  filter(type=="sample") %>% 
  select(ID) %>% 
  write_csv("../missing_samples_to_download.csv", col_names=FALSE)

# experiments
m_dat2 %>% 
  filter(type=="experiment") %>% 
  select(ID) %>% 
  write_csv("../missing_experiments_to_download.csv", col_names=FALSE)

