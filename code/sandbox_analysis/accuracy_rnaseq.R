require('tidyverse')
MIN.READS <- 100000
prefix <- "rat"
metadata_sex <- read.csv(sprintf("data/01_metadata/%s_rnaseq_metadata_sex.csv", prefix),
                         stringsAsFactors = FALSE) %>% 
  as_tibble() %>% 
  select(acc, mapped_sex) %>% 
  unique() 
table(metadata_sex$mapped_sex)
nrow(metadata_sex) 

exp_sample_nofilt <- read_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample.csv", prefix))
# mouse: 433339 samples, 6651 studies
exp_sample_nofilt2 <- exp_sample_nofilt %>% 
  semi_join(metadata_sex, by=c("sample_acc"="acc")) 
# mouse: 359142 samples, 6651 studies

exp_study1 <- exp_sample_nofilt2 %>% left_join(metadata_sex, by=c("sample_acc"="acc"))
exp_study2 <- exp_study1 %>% 
  group_by(study_acc) %>%
  summarize(num_samples=n(), 
            num_f=sum(mapped_sex=="female"),
            num_m=sum(mapped_sex=="male"))
exp_study3 <- exp_study2 %>% mutate(
  study_type=case_when(
    (num_f==0 & num_m==0) ~ "unknown",
    (num_f==0 & num_m != 0) ~ "male-only",
    (num_f!=0 & num_m == 0) ~ "female-only",
    TRUE ~ "mixed"
  )
)

exp_sample <- read_csv(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix))
exp_study <- exp_sample %>% 
  left_join(metadata_sex, by=c("sample_acc"="acc")) %>%
  mutate(qc_pass=(present & (num_reads >= MIN.READS))) 
# 121192 out of 359142, 3808 do not pass QC
# 6463 studyu



table(exp_study$qc_pass)

sex_lab <- metadata_sex %>% 
  filter(mapped_sex %in% c("male", "female")) 
  
sex_lab2 <- exp_sample %>% 
  left_join(sex_lab %>% rename(sample_acc=acc)) %>%
  
study_samp2 <- exp_study %>% 
  select(study_acc, sample_acc, mapped_sex, qc_pass) %>%
  rename(sex=mapped_sex)

study_counts_sex <- study_samp2 %>%
  group_by(study_acc) %>% 
  summarize( num_f=sum(sex=="female", na.rm=TRUE), 
             num_m=sum(sex=="male", na.rm=TRUE), num_samples=n(), 
             num_qc=sum(qc_pass)) %>%
  mutate(study_sex=case_when(
    num_samples < 10 ~ "small", # less than 10 samples
    num_f==num_samples ~ "female-only",
    num_m==num_samples ~ "male-only", 
    (num_m > 5 & num_f > 5) ~ "mixed sex",
    (num_f + num_m == 0) ~ "unknown",
    TRUE ~ "other")) # not all of the samples pass

table(study_counts_sex$study_sex)
mixed_sex_studies <- study_counts_sex %>% filter(study_sex=="mixed sex")
f_studies <- study_counts_sex %>% filter(study_sex=="female-only")
m_studies <- study_counts_sex %>% filter(study_sex=="male-only")

# load sex labels
imputed_sl <- read_csv(sprintf("data/02_labeled_data/%s_rnaseq_sl.csv", prefix)) %>%
  mutate(imputed_sex=ifelse(pred > 0.5, "male", "female"))

combined_lab <- left_join(study_samp2, imputed_sl, by=c("sample_acc"="id"))
combined_lab2 <- combined_lab %>% filter(!is.na(sex) & sex %in% c("female", "male") & qc_pass==TRUE) %>%
  select(sample_acc, sex, imputed_sex) %>% unique()
table(is.na(combined_lab2$imputed_sex)) 
# there are 39 w missing imputed labels, but 42 are missing in mouse
sum(combined_lab2$imputed_sex==combined_lab2$sex, na.rm=TRUE)/sum(!is.na(combined_lab2$imputed_sex))
# 90.0% concordance in human, 93.1% concordance in mouse 

# concordance in single-sex
f_samp <- combined_lab %>% filter(qc_pass==TRUE & !is.na(sex)) %>%
  filter(study_acc %in% f_studies$study_acc) %>%
  select(sample_acc, sex, imputed_sex) %>%
  unique()
table(is.na(f_samp$imputed_sex)) 

sum(f_samp$sex==f_samp$imputed_sex, na.rm=TRUE)/nrow(f_samp)

m_samp <- combined_lab %>% filter(qc_pass==TRUE & !is.na(sex) & sex=="male") %>%
  filter(study_acc %in% m_studies$study_acc) %>%
  select(sample_acc, sex, imputed_sex) %>%
  unique()
table(is.na(m_samp$imputed_sex)) 

sum(m_samp$sex==m_samp$imputed_sex, na.rm=TRUE)/nrow(m_samp)


# concordance in mixed sex HC
mm <- combined_lab %>% filter(qc_pass==TRUE & !is.na(sex) & sex %in% c("female", "male")) %>%
  filter(study_acc %in% mixed_sex_studies$study_acc) %>%
  select(study_acc, sample_acc, sex, imputed_sex) %>%
  unique()

