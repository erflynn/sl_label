
require('tidyverse')
require('data.table')
options(stringsAsFactors = FALSE)
# set up data for compendia train/test for SL, CL

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
data_type <- args[2]

if (data_type=="microarray"){
  mapping <- fread(sprintf("data/01_metadata/%s_compendia_dedup_mapping.csv", prefix, data_type), data.table = FALSE) %>%
    rename(study_acc=study_str)
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


set.seed(114)

no_cell <- metadata %>% filter(cl_line=="" | cl_line=="--")  %>% select(acc)

cell <- metadata %>% filter(cl_line!="" & cl_line!="--") %>% 
  right_join(refine_mapped %>% select(gsm, accession), by=c("acc"="gsm")) %>%
  filter(!is.na(cl_line)) %>% select(acc, cl_line, accession)

pos_cl <- cell %>% left_join(mapping, by=c("acc"="sample_acc")) %>% select(-cl_line, -accession) 
neg_cl <- no_cell %>% left_join(mapping, by=c("acc"="sample_acc"))  

refine_mapped %>% 
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

if (data_type=="microarray"){
  cl_annot_all <- inner_join(both_cell, study_cl2, by="gsm") %>%
    inner_join(metadata %>% select(acc, f_idx, idx), by=c("gsm"="acc")) %>%
    select(gse, gsm, is_sample_cl, is_study_cl, f_idx, idx)
} else {
  cl_annot_all <- inner_join(both_cell, study_cl2, by="gsm") %>%
    inner_join(metadata %>% select(acc), by=c("gsm"="acc")) %>%
    select(gse, gsm, is_sample_cl, is_study_cl)
}


cl_annot_all %>% write_csv(sprintf("data/02_labeled_data/%s_%s_cl_annot_info.csv", prefix, data_type2))
neg_cl <- cl_annot_all %>% filter(!is_sample_cl & !is_study_cl) %>% unique() #400k
pos_cl <- cl_annot_all %>% filter(is_sample_cl) %>% unique() #88k

pos_cl_drop <- drop_greater_n(pos_cl %>% select(acc, study_acc) %>% rename(grp=study_acc), 5, 2000) 
neg_cl_drop <- drop_greater_n(neg_cl %>% select(acc, study_acc) %>% rename(grp=study_acc), 5, 2000)

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
  pos_cl_drop <- pos_cl_drop %>% rename(study_acc=grp, sample_acc=acc)
  neg_cl_drop <- neg_cl_drop %>% rename(study_acc=grp, sample_acc=acc)
  
}

write_csv(neg_cl_drop, sprintf("data/03_train/%s_%s_cl_neg_sample.csv", prefix, data_type))
write_csv(pos_cl_drop, sprintf("data/03_train/%s_%s_cl_pos_sample.csv", prefix, data_type))



# ----- grab sex labeling info ----- #
sl <- read.csv(sprintf("data/01_metadata/%s_%s_metadata_sex.csv", prefix, data_type)) %>% unique()
sl2 <- sl %>% filter(sex==mapped_sex & sex %in% c("male", "female")) %>%
  rename(sample_acc=acc) %>%
  mutate(class=case_when(sex=="female" ~ 0, sex=="male" ~ 1)) %>% 
  select(-sex, -mapped_sex) 

no_cl <- metadata %>% filter(cl_line=="" | cl_line=="--")

map_sl <- sl2 %>% left_join(mapping) %>% unique() %>% rename(grp=study_acc)

map_sl2 <- map_sl %>% semi_join(no_cl, by=c(sample_acc="acc"))
table(map_sl2$class)

set.seed(114)
pos_drop <- drop_greater_n(map_sl %>% filter(class==1), 5, 2000) 
neg_drop <- drop_greater_n(map_sl %>% filter(class==0), 5, 2000)

# make output happy for next step
if (data_type=="microarray"){
  pos_drop <- pos_drop %>% rename(acc=sample_acc) %>% select(-class, -grp) %>% left_join(metadata %>% select(acc, f_idx, idx)) %>% arrange(idx)
  neg_drop <- neg_drop %>% rename(acc=sample_acc) %>% select(-class, -grp) %>% left_join(metadata %>% select(acc, f_idx, idx)) %>% arrange(idx)
} else {
  pos_drop <- pos_drop %>% rename(study_acc=grp) %>% select(study_acc, sample_acc)
  neg_drop <- neg_drop %>% rename(study_acc=grp) %>% select(study_acc, sample_acc)
  
}

write_csv(neg_drop, sprintf("data/03_train/%s_%s_sex_neg_sample.csv", prefix, data_type))
write_csv(pos_drop, sprintf("data/03_train/%s_%s_sex_pos_sample.csv", prefix, data_type))
