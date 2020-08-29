require('data.table')
require('tidyverse')
options(stringsAsFactors = FALSE)
prefix <- "human"
data_type <- "compendia"
MIN.READS <- 100000
# read in all the appropriate data

if (data_type=="compendia"){
  metadata <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", prefix))
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample.csv",prefix), data.table = FALSE)
} else {
  metadata <- read.csv(sprintf("data/01_metadata/%s_%s_sample_metadata.csv", prefix, data_type))
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix), data.table = FALSE)
  mapping <- mapping %>% filter(present & num_reads >= MIN.READS)
  
}

no_cell <- metadata %>% filter(cl_line=="")
refine_mapped <- read_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl.csv", 
                                  prefix, data_type))
refine_mapped2 <- read_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl_part2.csv",
                                   prefix, data_type))

cl_stud <- fread(sprintf("data/02_labeled_data/%s_%s_study_cl.csv", prefix, data_type), data.table = FALSE)


cell <- refine_mapped %>% 
  select(gsm, accession) %>% 
  bind_rows(refine_mapped2 %>% select(gsm, accession))

both_cell <- no_cell %>% anti_join(cell, by=c("acc"="gsm")) %>% 
  select(acc) %>%  mutate(cl="", accession="") %>%
  rename(gsm=acc) %>% select(accession, cl, gsm) %>%
  bind_rows(cell) %>%
  rename(cl_descript_sample=cl, cl_acc_sample=accession) %>%
  mutate(is_sample_cl=ifelse(cl_acc_sample=="" & cl_descript_sample=="", FALSE, TRUE))

# we also want the metadata to -NOT- have CL mentions
study_sample_cl <- full_join(cl_stud %>% select(gse, accession), mapping, by=c("gse"="study_acc"))
study_cl2 <- study_sample_cl %>% select(gse, sample_acc, everything()) %>% 
  rename(cl_acc_study=accession, gsm=sample_acc) %>%
  mutate(is_study_cl=ifelse(is.na(cl_acc_study), FALSE, TRUE))

# randomly sample some of these

cl_annot_all <- inner_join(both_cell, study_cl2, by="gsm") %>% 
  inner_join(metadata %>% select(acc, f_idx, idx), by=c("gsm"="acc")) %>%
  select(gse, gsm, is_sample_cl, is_study_cl, f_idx, idx) 

cl_annot_all %>% write_csv("data/02_labeled_data/human_cl_annot_info.csv")

neg <- cl_annot_all %>% filter(!is_sample_cl & !is_study_cl) %>% unique() #400k
pos <- cl_annot_all %>% filter(is_sample_cl) %>% unique() #88k

# numbers by study 
pos %>% group_by(gse) %>% count() %>% arrange(desc(n))
neg %>% group_by(gse) %>% count() %>% arrange(desc(n))

# select max 5 per study
study_filt <- function(df){
  df1 <- df %>% group_by(gse) %>% filter(n() >= 5) %>% sample_n(5)
  df2 <- df %>% group_by(gse) %>% filter(n() < 5) 
  
  # remove duplicates
  df3 <- rbind(df1, df2) %>% ungroup(gse) %>% select(-gse) %>%
  return()
}

pos2 <- pos %>% study_filt()
neg2 <- neg %>% study_filt()


# there are *tons* of samples...
# now group by study

# training/testing

set.seed(304)
neg_ex <- sample(1:nrow(neg2), 6000, replace=FALSE)
pos_ex <- sample(1:nrow(pos2), 6000, replace=FALSE)
not_cell_train <- neg2[neg_ex,]
cell_train <- pos2[pos_ex,]



pos_samp <- cell_train %>% rename(acc=gsm) %>% select(acc, f_idx, idx) %>% 
  arrange(idx) %>% unique() 
pos_samp %>% write_csv("data/03_train/human_cell_pos_sample.csv")
neg_samp <- not_cell_train %>% rename(acc=gsm) %>% select(acc, f_idx, idx) %>% 
  arrange(idx) %>% unique() 

neg_samp %>% write_csv("data/03_train/human_cell_neg_sample.csv")
