#
#
# 1. what are the missing files?
# - resubmit and rerun download for these
# 2. get a complete triple dataset
# 3. descriptions

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




# ---- grab run to sample mapping ----- #
my_runs = list.files(pattern="run_info*")
run_dat <- do.call(rbind, lapply(my_runs, 
                                 read_tsv, col_type="ccccc", col_names=FALSE))
run_dat2 <- run_dat %>% as_tibble()
colnames(run_dat2) <- c("run", "descript", "runs", "samples", "experiments")
run_mapping <- run_dat2 %>% select(run, samples) 

# one multi-maps --> perhaps skip?
run_mapping %>% filter(str_detect(samples, ";"))

# ----- load sample data ------ #
my_samples = list.files(pattern="sample_attr*")
sample_dat <- do.call(rbind, lapply(my_samples, 
                                 read_tsv, col_type="ccc", col_names=FALSE))
sample_dat1 <- sample_dat %>% as_tibble() 
colnames(sample_dat1) <- c("sample", "key", "value")
sample_dat2 <- sample_dat1 %>%
  mutate(across(key:value, tolower)) 

sg_sample_attr <- sample_dat2 %>%
  filter(str_detect(key, "sex|gender")) 

cl_sample_attr <- sample_dat2 %>%
  filter(str_detect(key, "cell")) 

trt_sample_attr <- sample_dat2 %>%
  filter(str_detect(key, "treatment")) 


run_to_attr <- run_mapping %>% inner_join(sg_sample_attr, by=c("samples"="sample"))
run_to_attr %>% write_csv("data/rnaseq_metadata_sex_sra_run_to_attr.csv")

run_to_attr %>% group_by(key) %>% count() %>% arrange(desc(n))
run_to_attr %>% group_by(value) %>% count() %>% arrange(desc(n))

run_to_attr %>% 
  distinct(run, value) %>%
  group_by(value) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  filter(value %in% c("male", "female", "m", "f")) %>%
  ungroup() %>%
  summarize(sum(n))

read_csv("data/01_metadata/combined_human_mouse_meta.csv") %>%
  filter(data_type=="rnaseq") %>%
  group_by(metadata_sex) %>%
  count() 