# grabs the metadata
# removes rnaseq samples from the compendia
# performs filtering, checks
# constructs the "list_samples.csv"

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

sample_metadata %>%
  write_csv("data/01_sample_lists/list_samples.csv")

# put together study metadata
study_metadata <- base_df %>%
  pmap(function(organism, data_type){
    loadStudyMetadata(organism, data_type)}) %>%
  bind_rows() %>% 
  addSrcCol(study_acc)

study_metadata %>%
  write_csv("data/01_sample_lists/list_studies.csv")


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

sample_to_study %>%
  write_csv("data/01_sample_lists/sample_to_study.csv")

# ------------------------------------------------------- #

# //TODO double check assumptions abt overlaps in coverage



# # put it all together in one file
# human_compendia <- loadSampleMetadata("human", "microarray")
# mouse_compendia <- loadSampleMetadata("mouse", "microarray")
# human_rnaseq <- loadSampleMetadata("human", "rnaseq")
# mouse_rnaseq <- loadSampleMetadata("mouse", "rnaseq")
# 
# 
# # remove the RNA-seq samples from the microarray files
# human_rnaseq_in_compendia <- human_compendia %>% 
#   filter(str_detect(sample_acc, "^[E|D|S]RR"))
# mouse_rnaseq_in_compendia <- mouse_compendia %>% 
#   filter(str_detect(sample_acc, "^[E|D|S]RR"))
# print(nrow(human_rnaseq_in_compendia))
# print(nrow(mouse_rnaseq_in_compendia))
# 
# human_microarray <- human_compendia %>%
#   filter(!str_detect(sample_acc, "^[E|D|S]RR"))
# mouse_microarray <- mouse_compendia %>%
#   filter(!str_detect(sample_acc, "^[E|D|S]RR"))
# print(nrow(human_microarray))
# print(nrow(mouse_microarray))
# 
# print(nrow(human_rnaseq))
# print(nrow(mouse_rnaseq))
# 
# # make sure none of the rnaseq samples in compendia are "new"
# stopifnot(length(setdiff(human_rnaseq_in_compendia$sample_acc, human_rnaseq$sample_acc))==0)
# stopifnot(length(setdiff(mouse_rnaseq_in_compendia$sample_acc, mouse_rnaseq$sample_acc))==0)
# 
# # put everything together
# sample_list <- map_df(list(human_microarray, mouse_microarray, 
#                            human_rnaseq, mouse_rnaseq), bind_rows)
# 
# 
# 
# # ---- grab the study metadata ---- #
# 
# 
# # ---- put together study list ---- #
# # study, organism, data_type, resource
# 
# 
# 
# 
# 
# 
# 
# study_metadata 
# 
# 
# # list of rnaseq studies
# # list of GEO
# # list of AE
# 
# 
# # --- put together exp-to-sample list --- #
# # study, sample, organism, data_type
# 



