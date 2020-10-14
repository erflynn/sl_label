# 02_get_sample_list.R
# E Flynn
# Last updated: 10/13/2020
#
# grabs the refine-bio metadata
# removes rnaseq samples from the compendia
# 
# constructs the "list_samples.csv", "list_studies.csv"

require('data.table')
require('tidyverse')

# ---- helpful functions ---- #
loadSampleMetadata <- function(my_organism, my_data_type){
  sample_metadata <- fread(sprintf("data/01_sample_lists/rb_metadata/%s_%s_sample_metadata.csv", 
                                   my_organism, my_data_type),
                           data.table=FALSE, stringsAsFactors=FALSE) %>% 
    as_tibble() %>%
    rename(sample_acc=acc) %>%
    select(sample_acc, platform) %>%
    mutate(organism=my_organism,
           data_type=my_data_type)
  
  # make sure there is one row per sample
  stopifnot(length(unique(sample_metadata$sample_acc))==nrow(sample_metadata))

  return(sample_metadata %>% select(sample_acc, organism, data_type, platform)) 
}


loadStudyMetadata <- function(my_organism, my_data_type){
  study_metadata <- fread(sprintf("data/01_sample_lists/rb_metadata/%s_%s_experiment_metadata.csv", 
                                  my_organism, my_data_type),
                          data.table=FALSE, stringsAsFactors=FALSE) %>% 
    as_tibble() %>%
    select(study_acc) %>%
    mutate(organism=my_organism,
           data_type=my_data_type)
  
  # make sure there is one row per study
  stopifnot(length(unique(study_metadata$study_acc))==nrow(study_metadata))
  
  return(study_metadata) 
}


addSrcCol <- function(df, acc) {
  df2 <- df %>% mutate("src"=case_when(
    str_detect({{acc}}, "^[E|S|D]R") ~ "SRA",
    str_detect({{acc}}, "^GS") ~ "GEO",
    TRUE ~ "ArrayExpress"
  ))
  # stop if any RNA-seq samples are GEO or arrayExpress
  stopifnot(df2 %>% 
              filter(data_type=="rnaseq" & 
                       src %in% c("GEO", "ArrayExpress")) %>%
              nrow()==0)
  return(df2)
}  

# ----- load the data ----- #

base_df <- list(
  organism=c("human", "mouse"), 
  data_type=c("microarray","rnaseq"))%>%
  cross_df()

# put together sample metadata 
sample_metadata <-base_df %>%
  pmap(function(organism, data_type){
    loadSampleMetadata(organism, data_type)}) %>%
  bind_rows() %>% 
  addSrcCol(sample_acc)

sample_metadata2 <- sample_metadata %>%
  filter(src != "SRA" | (src=="SRA" & data_type=="rnaseq"))
stopifnot(nrow(sample_metadata2)==length(unique(sample_metadata2$sample_acc)))


# put together study metadata
study_metadata <- base_df %>%
  pmap(function(organism, data_type){
    loadStudyMetadata(organism, data_type)}) %>%
  bind_rows() %>% 
  addSrcCol(study_acc)

study_metadata2 <- study_metadata %>%
  filter(src != "SRA" | (src=="SRA" & data_type=="rnaseq"))

# unique id = study/organism 
study_metadata2 %>% distinct(study_acc, organism) %>% nrow() == nrow(study_metadata)



# put together sample/study mapping
sample_to_study <- base_df %>%
  pmap(function(organism, data_type){
    read_csv(sprintf("data/01_sample_lists/rb_metadata/%s_%s_exp_to_sample.csv", 
                     organism, data_type), col_types="cc") %>%
      mutate(organism=organism,
             data_type=data_type)
  }) %>%
  bind_rows() %>% 
  addSrcCol(sample_acc)

sample_to_study2 <- sample_to_study %>%
  filter(src != "SRA" | (src=="SRA" & data_type=="rnaseq")) 

# row id: sample, study, organism
stopifnot(sample_to_study2 %>% distinct(sample_acc, study_acc, organism) %>% nrow()==nrow(sample_to_study2))
# ------------------------------------------------------- #

# assumptions:
# 1. all studies are in samp_to_study
stopifnot(length(setdiff(study_metadata2$study_acc, sample_to_study$study_acc))==0)

# 2. all samples are in samp_to_study
stopifnot(length(setdiff(sample_metadata2$sample_acc, sample_to_study$sample_acc))==0)

# 3. all SRA samples in the compendia are also in RNA-seq
rnaseq_in_compendia <- sample_metadata %>% filter(src == "SRA" & data_type=="microarray")
rnaseq_samples <- sample_metadata %>% filter(data_type=="rnaseq")
stopifnot(length(setdiff(rnaseq_in_compendia$sample_acc, rnaseq_samples$sample_acc))==0)

# -------------- #
# write it out
sample_metadata2 %>%
  write_csv("data/01_sample_lists/list_samples.csv")

study_metadata2 %>%
  write_csv("data/01_sample_lists/list_studies.csv")

sample_to_study2 %>%
  write_csv("data/01_sample_lists/sample_to_study.csv")

