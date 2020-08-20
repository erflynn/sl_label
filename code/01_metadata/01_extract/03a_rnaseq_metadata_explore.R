
require('tidyverse')
my_runs = list.files(pattern="run_info*")
run_dat <- do.call(rbind, lapply(my_runs, 
                                 read_tsv, col_type="ccccc", col_names=FALSE))
run_dat2 <- run_dat %>% as_tibble()
colnames(run_dat2) <- c("run", "descript", "runs", "samples", "experiments")
run_mapping <- run_dat2 %>% select(run, samples) 

# one multi-maps --> perhaps skip?
run_mapping %>% filter(str_detect(samples, ";"))


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