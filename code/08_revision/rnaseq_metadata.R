
library(tidyverse)
rnaseq_h <- read_csv("data/01_sample_lists/rb_metadata/human_rnaseq_experiment_metadata.csv")
single_cell_studies <- rnaseq_h %>% 
  mutate(across(c(title, description), tolower)) %>%
  filter(str_detect(title, "single cell|scrna"),
         str_detect(description, "single cell|scrna")) # 153

rnaseq_m <- read_csv("data/01_sample_lists/rb_metadata/mouse_rnaseq_experiment_metadata.csv")
single_cell_studies_m <- rnaseq_m %>% 
  mutate(across(c(title, description), tolower)) %>%
  filter(str_detect(title, "single cell|scrna"),
         str_detect(description, "single cell|scrna")) # 187


comb_metadata <- read_csv("data/sample_metadata_filt.csv", col_types="cccccdldcc")
to_remove <- comb_metadata %>% 
  filter(study_acc %in% c(single_cell_studies_m$study_acc, single_cell_studies$study_acc)) %>%
  filter(!is.na(present))

to_keep <- comb_metadata %>% filter(data_type=="rnaseq" & !is.na(present)) %>%
  filter(!study_acc %in% c(single_cell_studies_m$study_acc, single_cell_studies$study_acc)) 

plot(density(to_keep$num_reads))
plot(density(to_remove$num_reads))

rnaseq_ds <- to_keep %>% mutate(keep=TRUE) %>% bind_rows(to_remove %>% mutate(keep=FALSE)) %>%
  select(-platform, -present, -data_type, -study_acc) %>%
  pivot_wider(names_from="label_type", values_from="sex_lab") 
ggplot(rnaseq_ds, aes(x=num_reads, fill=expression))+
  geom_histogram(alpha=0.9, bins=100)+
  facet_grid(organism~keep, scales="free")+
  theme_bw()+
  xlim(c(0, 40000000))+
  xlab("number of reads")+
  ylab("number of samples")

# 1/3 of these are unlabeled
to_keep %>% 
  filter(label_type=="expression" & num_reads > 5000000) %>% 
  group_by(sex_lab=="unlabeled") %>% 
  count() %>%
  ungroup() %>%
  mutate(frac=n/sum(n))

# 2/3 of these are unlabeled
to_remove %>% 
  filter(label_type=="expression") %>% 
  group_by(sex_lab=="unlabeled") %>% 
  count() %>%
  ungroup() %>%
  mutate(frac=n/sum(n))


