# Process our metadata sex labels and then compare to MetaSRA's
#
# TODO:
# - what is the reason it doesn't overlap

require('tidyverse')
# read in our sex labels
sg_tab2 <- read_csv("data/01_metadata/mapped_sl_all.csv")

run_to_sample <- read_csv("data/sra_run_to_sample.csv")

# read in the combined metadata 
comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", 
                          col_types="cccccccdcd")
# collapse
sg_tab3 <- sg_tab2 %>% 
  filter(mapped_sex!="") %>%
  group_by(sample) %>%
  summarise(mapped_sex=paste(unique(mapped_sex), collapse=";"),
            type_key=paste(unique(type_key), collapse=";")) 
sg_tab4 <- sg_tab3 %>%
  mutate(mapped_sex=ifelse(str_detect(mapped_sex, ";"),"unknown", mapped_sex))
stopifnot((sg_tab4 %>% distinct(sample) %>% nrow())==nrow(sg_tab4))            

# add a run column
sg_run <- run_to_sample %>% inner_join(sg_tab4, by=c("samples"="sample")) %>%
  rename(sample=samples, sex=mapped_sex, key_type=type_key) 
stopifnot((sg_run %>% distinct(run) %>% nrow())==nrow(sg_run))            
sg_run %>% write_csv("data/01_metadata/run_mapped_rnaseq_sex.csv")


# get the list of all samples in our dataset
rnaseq_dat <- comb_metadata %>% filter(data_type=="rnaseq")
rnaseq2 <- rnaseq_dat %>% 
  select(sample_acc, organism, metadata_sex) %>% 
  left_join(sg_run, by=c("sample_acc"="run")) %>%
  mutate(sex=ifelse(is.na(sex), "unknown", sex))
table(rnaseq2$metadata_sex, rnaseq2$sex)

sum(rnaseq2$metadata_sex=="unknown")-sum(rnaseq2$sex=="unknown") # 45.7k
# better :)

# not organism-specific
rnaseq2 %>% 
  group_by(organism) %>% 
  summarise(ms=sum(metadata_sex=="unknown"), ss=sum(sex=="unknown"), n=n()) %>%
  mutate(diff=ms-ss) %>%
  mutate(across(c(ss, ms, diff), ~./n))

#sum(rnaseq2$sex=="unknown")/nrow(rnaseq2) # 0.7600228
#sum(rnaseq2$metadata_sex=="unknown")/nrow(rnaseq2) # 0.8376703

rnaseq2 %>% 
  filter(metadata_sex=="unknown" & sex!="unknown") %>%
  common_col(key_type) # 44482 come from std, 1373 come from rescue

rnaseq2 %>% 
  filter(metadata_sex=="unknown" & sex!="unknown") %>%
  sample_n(10) %>%
  left_join(sg_tab2, by="sample") %>% select(key, value, sex)

# --- Compare to MetaSRA --- #
load("data/10_comparison/metasra_data.RData") # --> mapped_df

metasra_sl <- mapped_df %>% 
  filter(term_id %in% c("UBERON:0003100", "UBERON:0003101")) %>%
  mutate(metasra_sex=case_when(term_id=="UBERON:0003100" ~ "female",
                               term_id=="UBERON:0003101" ~ "male"))  %>%
  group_by(sample_accession) %>%
  summarise(metasra_sex=paste(unique(metasra_sex), collapse=";")) %>%
  mutate(metasra_sex=ifelse(str_detect(metasra_sex, ";"), "mixed", metasra_sex)) %>%
  ungroup() 
nrow(metasra_sl) # 153417
nrow(sg_tab4) # 127755

# NOTE - metasra has no mouse! 
# but somehow has more data? but less for our samples

# put together the data for comparison
comparison_tab <- metasra_sl %>%   
  left_join(run_to_sample, by=c("sample_accession"="samples")) %>%
  right_join(rnaseq2, by=c("run"="sample_acc")) %>%
  mutate(metasra_sex=ifelse(is.na(metasra_sex), "unknown", metasra_sex)) %>%
  rename(rb_sex=metadata_sex, mapped_sex=sex, sample_acc=run) %>%
  select(-sample_accession, -sample,-key_type) %>%
  select(sample_acc, organism, everything())
comparison_tab %>% 
  filter(organism=="human") %>% 
  summarise(ms=sum(metasra_sex=="unknown"), ss=sum(mapped_sex=="unknown"), n=n()) %>%
  mutate(diff=ms-ss) %>%
  mutate(across(c(ss, ms, diff), ~./n))

# create a confusion matrix of MetaSRA (rows) vs our labels (cols)
conf_mat_human <- comparison_tab %>% 
  filter(organism=="human") %>% 
  group_by(metasra_sex, mapped_sex) %>% 
  count() %>% 
  ungroup() %>%
  pivot_wider(names_from="mapped_sex", values_from="n", values_fill=0) %>%
  select(metasra_sex, female, male, mixed, unknown)
conf_mat_human %>% write_csv("tables/metasra_conf_mat_sex.csv")
