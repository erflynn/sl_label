#
# Pull the RNA-seq sample attributes after parsing XML from ENA

require('tidyverse')

# ---- grab run to sample mapping ----- #
my_runs = list.files(path="data/sra_out", pattern="run_info*")
run_dat <- do.call(rbind, lapply(my_runs, 
                                 read_tsv(sprintf("data/sra_out/%s", x), 
                                          col_type="ccccc", 
                                          col_names=FALSE)))
run_dat2 <- run_dat %>% as_tibble()
colnames(run_dat2) <- c("run", "descript", "runs", "samples", "experiments")
run_mapping <- run_dat2 %>% select(run, samples) 
run_mapping %>% write_csv("data/sra_run_to_sample.csv")


# --------- grab sample attrib --------- #
my_samples = list.files(path="data/sra_out", pattern="sample_attr*")
sample_dat <- do.call(rbind, lapply(my_samples, function(x)
  read_tsv(sprintf("data/sra_out/%s", x), col_type="ccc", 
           col_names=FALSE)))

# cleaning
sample_dat1 <- sample_dat %>% as_tibble() 
colnames(sample_dat1) <- c("sample", "key", "value")
sample_dat2 <- sample_dat1 %>%
  mutate(across(key:value, tolower)) 

# remove ENA data for now
sample_dat3 <- sample_dat2 %>% 
  filter(!str_detect(key, "ena-"))
rm(sample_dat, sample_dat1, sample_dat2)
save(sample_dat3, file="data/sample_to_attr_sm.RData")