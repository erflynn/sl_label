# Check the metadata sex labels for RNA-seq data
#
# There are fewer RNA-seq sample sex labels than microarray.
# This doesn't quite make sense --
# and it looks like the issue is metadata parsing! :/
# TODO: redownload and redo all of this


require('tidyverse')
load("data/dates/sample_attr_sra.RData") # -- sample_ds with attributes
load("data/dates/run_sample_sra.RData") # -- h_runs_to_sample, m_runs_to_sample

ds <- sample_ds %>%
  filter(!is.na(sample_attribute)) %>%
  as_tibble() %>%
  separate_rows(sample_attribute, sep=" \\|\\| ") %>%
  separate(sample_attribute, into=c("key", "value"), sep=": ", extra="merge")

unique(ds$key)
ds %>% group_by(key) %>% count() %>% arrange(desc(n))
keys <- (ds %>% group_by(key) %>% count() %>% arrange(desc(n)))$key

sex_lab_keys <- keys[str_detect(sapply(keys,tolower), "sex|gender")]

sex_lab_sample_attr <- ds %>% 
  filter(key %in% sex_lab_keys) %>%
  mutate(value=tolower(value)) 

cleaned_sl <- sex_lab_sample_attr %>% 
  select(value) %>% 
  unique() %>%
  group_by(value) %>%
  mutate(sex_lab=case_when(
    is.na(value) ~ "unknown",
    value %in% c("not applicable","not determined", "not collected","undetermine" ,
                 "not available", "missing", "none", "asexual") ~ "unknown",
    value %in% c("", "na", "n/a", "-", "--", "no", "?", "u", "dk") ~ "unknown",
    str_detect(value, "mixed|both") ~ "mixed",
    any(str_split(value, ",| ")[[1]]=="female") & 
      any(str_split(value, ",| ")[[1]]=="male") ~ "mixed",
    any(str_split(value, ",| ")[[1]]=="female") ~ "female",
    any(str_split(value, ",| ")[[1]]=="male") ~ "male",
    value %in% c("b", "fmm", "mf", "fm", "ffm", "m and f") ~ "mixed",
    value %in% c("xy", "m", "mm") ~ "male",
    value %in% c("xx", "f") ~ "female",
    TRUE ~ value
  )) %>%
  ungroup()

sex_lab_sample_attr2 <- sex_lab_sample_attr %>% 
  left_join(cleaned_sl, by=c("value")) %>%
  rename(sample=sample_accession) %>%
  select(sample, sex_lab)

# mapping table
runs_to_sample_sl <- h_runs_to_sample %>% 
  bind_rows(m_runs_to_sample) %>%
  left_join(sex_lab_sample_attr2, by="sample") %>%
  rename(sample_acc=run, sra_sample_id=sample, sra_sample_sex=sex_lab) %>%
  unique() %>%
  group_by(sample_acc) %>%
  mutate(sra_sample_id=paste(sra_sample_id, collapse=";"),
         sra_sample_sex=paste(unique(sra_sample_sex), collapse=";")) %>%
  unique() %>%
  mutate(sra_sample_sex=ifelse(sra_sample_sex=="NA", "unknown", sra_sample_sex))

stopifnot(length(unique(runs_to_sample_sl$sample_acc))==nrow(runs_to_sample_sl))
table(runs_to_sample_sl$sra_sample_sex)



# read in what we have already!
comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")

rnaseq_w_sample <- comb_metadata %>% 
  filter(data_type=="rnaseq") %>%
  left_join(runs_to_sample_sl, by="sample_acc") 

table(is.na(rnaseq_w_sample$sra_sample_sex)) # 5110 missing
stopifnot(length(unique(rnaseq_w_sample$sample_acc))==nrow(rnaseq_w_sample))

rnaseq_w_sample2 <- rnaseq_w_sample %>%
  mutate(sra_sample_sex=ifelse(is.na(sra_sample_sex), "unknown", sra_sample_sex)) %>%
  mutate(sra_sample_id=ifelse(sra_sample_id=="NA", NA, sra_sample_id)) %>%
  select(-platform, -data_type) 
rnaseq_w_sample2 %>% group_by(metadata_sex, sra_sample_sex) %>% count()

# in SRA but not rb
rnaseq_w_sample2 %>% 
  filter(metadata_sex=="unknown" & sra_sample_sex != "unknown") %>%
  nrow() # 34522 (5.86%, 26.5% lab)

# in rb but not SRA
rnaseq_w_sample2 %>% 
  filter(metadata_sex!="unknown" & sra_sample_sex == "unknown") %>%
  nrow() # 56511 (9.60%, 43.4% lab)

rnaseq_w_sample2 %>% nrow() # 588931
rnaseq_w_sample2 %>% 
  filter(metadata_sex!= "unknown" |sra_sample_sex!="unknown") %>%
  nrow() # 130123 labeled in either

# matching count
rnaseq_w_sample2 %>% 
  filter(metadata_sex!= "unknown", sra_sample_sex==metadata_sex) %>%
  nrow() # 39090 (6.64% total, 30.0% labeled)

rnaseq_w_sample2 %>% 
  filter(metadata_sex!= "unknown",
         sra_sample_sex!="unknown", 
         sra_sample_sex!=metadata_sex) %>%
  nrow() # 0


# -- check the run attributes! -- #
# do these match the data?
load("data/dates/run_attr_sra.RData") # -- run_ds with attributes

run_attr <- run_ds %>%
  filter(!is.na(run_attribute)) %>%
  as_tibble() %>%
  separate_rows(run_attribute, sep=" \\|\\| ") %>%
  separate(run_attribute, into=c("key", "value"), sep=": ", extra="merge") %>%
  mutate(value=tolower(value))


run_keys <- (run_attr %>% group_by(key) %>% count() %>% arrange(desc(n)))$key
sex_lab_keys_run <- run_keys[str_detect(sapply(run_keys,tolower), "sex|gender")]

run_attr %>% 
  filter(key %in% sex_lab_keys_run) %>%
  group_by(value) %>% 
  count()

# there are v few of these!! how did refine-bio get these annotations?

# -- should I make sure the compendia sex labels are ok too? -- #