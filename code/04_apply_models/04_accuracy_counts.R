require('tidyverse')

prefix <- "human"
data_type <- "microarray"
ds <- "sex"

my_dat <- read_csv(sprintf("data/09_model_summary/%s_%s_%s_labels.csv", prefix, data_type, ds))
metadata_sl <- read.csv(sprintf("data/01_metadata/%s_%s_metadata_sex.csv", prefix, data_type)) %>% unique()

load(sprintf("data/07_model_dat/%s_%s_sex_train_mat.RData", prefix, data_type)) # --> X_train, X_test, Y_train, Y_test

if (data_type=="microarray"){
  mapping <- fread(sprintf("data/01_metadata/%s_compendia_dedup_mapping.csv", prefix, data_type), data.table = FALSE) %>%
    rename(study_acc=study_str)
  data_type2 <- "compendia"
  metadata <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", prefix))
} else {
  MIN.READS <- 100000
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix), data.table=FALSE)  %>%
    filter(present & num_reads >= MIN.READS) %>% select(-present)
  data_type2 <- data_type
  metadata <- read.csv(sprintf("data/01_metadata/%s_%s_sample_metadata.csv", prefix, data_type))
  
}



my_lab <- my_dat %>% mutate(ypred=round(pred)) %>% rename(sample_acc=id)

sl_comp <- metadata_sl %>% 
  filter(mapped_sex %in% c("male", "female")) %>% 
  select(-sex) %>% 
  rename(sample_acc=acc, sex=mapped_sex) %>%
  mutate(y=ifelse(sex=="female", 0,1)) %>%
  left_join(my_lab) %>% 
  mutate(in_train=(sample_acc %in% rownames(X_train))) %>%
  mutate(in_test=(sample_acc %in% rownames(X_test))) %>%
  filter(!is.na(ypred)) # why are so many missing?


get_frac <- function(ds){
  ds %>% group_by(y, ypred) %>% 
    count() %>% 
    ungroup() %>% 
    mutate(eq=(ypred==y)) %>%
    group_by(eq) %>% 
    summarize(n=sum(n)) %>% 
    ungroup() %>%
    mutate(tot=sum(n)) %>% 
    group_by(eq) %>%
    mutate(frac=n/tot)
}

# look at the amt not in the train/test
sl_comp %>% 
  filter(!in_train & !in_test) %>% 
  get_frac() # human RNA-seq: 91.6%, human microarray: 92.9%

# look at all
sl_comp %>% 
 get_frac() # 90.6%

# ok so in train isnt killing it
#  -- but what abt leakage from same studies?
# // TODO: what happens if we filter by study_str?




# // TODO: HC


metadata_sl2 <- mapping %>% left_join(metadata_sl, by=c("sample_acc"="acc"))

metadata_sl3 <- metadata_sl2 %>% group_by(study_acc) %>%
  summarize(num_f=sum(sex=="female"),
            num_m=sum(sex=="male"),
            num_tot=n())

f_only <- metadata_sl3 %>% filter(num_f>5 & num_f==num_tot)
m_only <- metadata_sl3 %>% filter(num_m>5 & num_m==num_tot)

f_dat <- sl_comp %>% semi_join(f_only %>% left_join(metadata_sl2) %>% select(sample_acc)) 
m_dat <- sl_comp %>% semi_join(m_only %>% left_join(metadata_sl2) %>% select(sample_acc)) 

f_dat %>% get_frac() # 99.2, 98.1
m_dat %>% get_frac() # 71.8, 91.1

f_dat %>% filter(!in_train & !in_test) %>% get_frac() # 99.4%
m_dat %>% filter(!in_train & !in_test) %>% get_frac() # 72.0%

# there are overall more m labeled as f than f labeled as m... hmm
table(sl_comp[,c("y", "ypred")])

# look at this by cell line?
sample_cl <- read_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl.csv", 
                                  prefix, data_type2))
sample_cl2 <- read_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl_part2.csv",
                                   prefix, data_type2))

study_cl <- fread(sprintf("data/02_labeled_data/%s_%s_study_cl.csv", prefix, data_type2), data.table = FALSE) %>% 
  left_join(metadata_sl2, by=c("gse"="study_acc"))



f_dat %>% filter(!sample_acc %in% sample_cl$gsm & 
                   !sample_acc %in% sample_cl2$gsm & 
                   !sample_acc %in% study_cl$sample_acc) %>% get_frac() # 99.2%
m_dat %>% filter(!sample_acc %in% sample_cl$gsm & 
                   !sample_acc %in% sample_cl2$gsm & 
                   !sample_acc %in% study_cl$sample_acc) %>% get_frac() 
# 83.0% -- so cell line is def part of it, but why is it still lower?

mismatch_m <- m_dat %>% filter(!sample_acc %in% sample_cl$gsm & 
                             !sample_acc %in% sample_cl2$gsm & 
                             !sample_acc %in% study_cl$sample_acc) %>%
   filter(y!=ypred) 

correct_m <- m_dat %>% filter(!sample_acc %in% sample_cl$gsm & 
                                 !sample_acc %in% sample_cl2$gsm & 
                                 !sample_acc %in% study_cl$sample_acc) %>%
  filter(y==ypred) 

# these look different!!
plot(density(mismatch_m$pred)) # this is more inward
plot(density(1-correct_m$pred))

m_dat %>% filter(pred < 0.2 | pred > 0.8) %>% get_frac() # 98.3% -- but this is only 1/3
f_dat %>% filter(pred < 0.2 | pred > 0.8) %>% get_frac() # 98.9% -- only 1/3

m_dat %>% filter(pred < 0.3 | pred > 0.7) %>% get_frac() # 97.2% -- ~ 1/2
f_dat %>% filter(pred < 0.3 | pred > 0.7) %>% get_frac() # 99.1% -- ~1/2


# // TODO: put all this together and summarize!


