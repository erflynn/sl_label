
# DEFINE CELL LINE VS NON-CELL LINE, then set up train/test

require('tidyverse')
require('data.table')
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
data_type <- args[2]

if (data_type=="microarray"){
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample.csv", prefix, data_type), data.table = FALSE) 
  data_type2 <- "compendia"
  metadata <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", prefix))
} else {
  MIN.READS <- 100000
  mapping <- fread(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix), data.table=FALSE)  %>%
    filter(present & num_reads >= MIN.READS)
  data_type2 <- data_type
  metadata <- read.csv(sprintf("data/01_metadata/%s_%s_sample_metadata.csv", prefix, data_type))
  
}


# drops samples such that no more than n samples are in a group
drop_greater_n <- function(ds, n_group, num_tot){
  ds2 <- ds %>% group_by(grp) %>% filter(n() >=n_group) %>% sample_n(n_group)
  ds3 <- ds %>% group_by(grp) %>% filter(n() < n_group) 
  ds4 <- rbind(ds2, ds3) %>% ungroup()
  if (nrow(ds4) <= num_tot){
    return(ds4)
  } else {
    ds5 <- ds4 %>% sample_n(num_tot)
    return(ds5)  
  }
}



# ----- grab cell labeling info ----- #
refine_mapped <- read_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl.csv", 
                                  prefix, data_type2))
refine_mapped2 <- read_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl_part2.csv",
                                   prefix, data_type2))

cl_stud <- fread(sprintf("data/02_labeled_data/%s_%s_study_cl.csv", prefix, data_type2), data.table = FALSE)

no_cell <- metadata %>% filter(cl_line=="" | cl_line=="--")  %>% select(acc)
some_cl <- metadata %>% filter(cl_line!="" | cl_line!="--")  

cl_lab_unique <- refine_mapped  %>% 
  bind_rows(refine_mapped2 %>% select(gsm, orig_str, accession)) %>% 
  unique()

xenograft_cl <- cl_lab_unique %>% filter(str_detect(orig_str, "xenograft"))
primary_cl <- cl_lab_unique %>% filter(str_detect(orig_str, "primary"))
stem_cl <- cl_lab_unique %>% filter(str_detect(orig_str, "stem|ipsc"))

meta_unique <- metadata %>% select(acc, cl_line, part, title) 

#unmapped_cl <- some_cl %>% anti_join(cl_lab_unique, by=c("acc"="gsm"))

mu_s <- meta_unique %>% left_join(cl_lab_unique, by=c("acc"="gsm")) %>% 
  group_by(acc)  %>% mutate(str=paste(c(cl_line, part, title), collapse=" ")) %>%
  ungroup() %>%
  mutate(str=tolower(str)) # SLOW
mu_s2 <- mu_s %>% mutate(source_type=case_when(
  str_detect(str, "xenograft") ~ "xenograft",
  str_detect(str, "stem cell|ipsc") ~ "stem_cell",
  !is.na(accession) & !str_detect(str, "primary") ~ "named_cl",
  !is.na(accession) & str_detect(str, "primary") ~ "other",
  str_detect(str, "tumor|cancer|carcinoma|melanoma|malignant") ~ "cancer",
  str_detect(str, "primary") ~ "primary_cells",
  str_detect(str, "cell line") ~ "unnamed_cl",
  str_detect(str, "cell|culture|passage")  ~ "other",
  is.na(cl_line) | cl_line %in% c("", " ", "--") ~ "tissue",
  TRUE ~ "other"
))

table(mu_s2$source_type)

# now add in the study_level
samp_w_study_cl <- study_cl %>% select(gse, accession) %>% unique() %>% left_join(mapping, by=c("gse"="study_acc")) 
table((mu_s2 %>% semi_join(samp_w_study_cl, by=c("acc"="sample_acc")))$source_type)
# ok most of these already assigned
mu_s3 <- mu_s2 %>% mutate(source_type=ifelse(source_type=="tissue" & acc %in% samp_w_study_cl$sample_acc, "other", source_type))
table((mu_s3 %>% semi_join(samp_w_study_cl, by=c("acc"="sample_acc")))$source_type)
sample_source_df <- mu_s3 %>% select(acc, source_type)
sample_source_df %>% write_csv(sprintf("data/02_labeled_data/%s_%s_sample_source.csv", prefix, data_type2))

### SAVE ###

set.seed(114)
### STOP ###

pos_cl <- sample_source_df %>% filter(source_type=="named_cl") %>% left_join(mapping, by=c("acc"="sample_acc")) %>% select(-source_type)
neg_cl <- sample_source_df %>% filter(source_type=="tissue")  %>% left_join(mapping, by=c("acc"="sample_acc")) %>% select(-source_type)
# note - these do overlap!!


# randomly sample some of these
pos_cl_drop <- drop_greater_n(pos_cl %>% select(acc, study_acc) %>% dplyr::rename(grp=study_acc), 5, 2000) 
neg_cl_drop <- drop_greater_n(neg_cl %>% select(acc, study_acc) %>% dplyr::rename(grp=study_acc), 5, 2000)

if (data_type=="microarray"){
  pos_cl_drop <- pos_cl_drop %>% 
    left_join(metadata %>% select(acc, f_idx, idx)) %>% 
    select(-grp) %>%
    unique() %>%
    arrange(idx)
  neg_cl_drop <- neg_cl_drop %>% 
    left_join(metadata %>% select(acc, f_idx, idx)) %>% 
    select(-grp) %>%
    unique() %>%
    arrange(idx)
} else {
  pos_cl_drop <- pos_cl_drop %>% dplyr::rename(study_acc=grp, sample_acc=acc)
  neg_cl_drop <- neg_cl_drop %>% dplyr::rename(study_acc=grp, sample_acc=acc)
}

pos_cl_drop %>% sample_n(20) %>% left_join(mu_s3)
neg_cl_drop %>% sample_n(20) %>% left_join(mu_s3)

write_csv(neg_cl_drop, sprintf("data/03_train/%s_%s_cl_neg_sample.csv", prefix, data_type))
write_csv(pos_cl_drop, sprintf("data/03_train/%s_%s_cl_pos_sample.csv", prefix, data_type))
