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
library('GEOmetadb')


# ---- helpful functions ---- #
loadSampleMetadata <- function(my_organism, my_data_type){
  sample_metadata <- read_csv(sprintf("data/01_sample_lists/rb_metadata/%s_%s_sample_metadata.csv", 
                                   my_organism, my_data_type),
                           col_types="cccccccccccc") %>% 
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
  study_metadata <- read_csv(sprintf("data/01_sample_lists/rb_metadata/%s_%s_experiment_metadata.csv", 
                                  my_organism, my_data_type), col_types="cccc") %>% 
    as_tibble() %>%
    #select(study_acc, date) %>%
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
  addSrcCol(sample_acc) %>%
  mutate(sample_acc=str_trim(sample_acc))


# put together study metadata
study_metadata <- base_df %>%
  pmap(function(organism, data_type){
    loadStudyMetadata(organism, data_type)}) %>%
  bind_rows() %>% 
  addSrcCol(study_acc) %>%
  distinct() 


# put together sample/study mapping
sample_to_study <- base_df %>%
  pmap(function(organism, data_type){
    read_csv(sprintf("data/01_sample_lists/rb_metadata/%s_%s_exp_to_sample.csv", 
                     organism, data_type), col_types="cc") %>%
      mutate(organism=organism,
             data_type=data_type)
  }) %>%
  bind_rows() %>% 
  filter(!is.na(sample_acc))  %>%
  mutate(sample_acc=str_trim(sample_acc))
# ------------------------------------------------------- #

# NEXT: label potential scRNA-seq data studies
study_metadata2 <- study_metadata %>% 
  mutate(across(c(title, description), tolower)) %>%
  mutate(scrna=(str_detect(title, "single cell|scrna") |
         str_detect(description, "single cell|scrna"))) # --> 1281 studies

scrna_likely <- sample_to_study %>% 
  semi_join(study_metadata2 %>% filter(scrna), by="study_acc") %>% 
  distinct(sample_acc) # 405,775 scrnaseq samples
sample_metadata2 <- sample_metadata %>%
  mutate(scrna=(sample_acc %in% scrna_likely$sample_acc))


# ----- CLEAN UP PLATFORM DATA ---- #


geo_samples <- unique(sample_metadata2 %>% filter(src=="GEO") %>% pull(sample_acc))
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
gsm_to_gpl <- dbGetQuery(con, sprintf("SELECT gsm, gpl FROM gsm WHERE gsm IN ('%s');", 
                                      paste(geo_samples, collapse="','")))
gpl_info <- dbGetQuery(con, sprintf("SELECT gpl, title FROM gpl WHERE gpl IN ('%s');", 
                                    paste(unique(gsm_to_gpl$gpl), collapse="','")))

# sample_metadata2 %>% 
#   filter(str_detect(platform, "zebrafish|bovine|rat|drosophila")) %>% 
#   left_join(gsm_to_gpl, by=c("sample_acc"="gsm")) %>%
#   left_join(gpl_info)

gpl_info2 <- gpl_info %>% 
  mutate(short_title=str_extract(title, "^\\[[A-z0-9\\-\\_ \\.]+\\]"))

gpl_info3 <- gpl_info2 %>% 
  group_by(title) %>%
  mutate(title2=ifelse(is.na(short_title),
    str_split(title, "\\[|\\(")[[1]][[1]], 
    str_split(title, "\\[|\\]|\\(")[[1]][[3]])) %>%
  mutate(title2=str_trim(str_replace(title2, "\\[|\\]|\\(", "")),
         short_title=str_replace(str_split(short_title, "\\]")[[1]][[1]], "\\[", "")) %>%
  ungroup() 

gpl_info4 <- gpl_info3 %>%
  group_by(title2) %>%
  mutate(title3= ifelse(
    str_detect(title2, "^-"), str_split(title, "- ")[[1]][[2]], title2)) %>%
  mutate(title3=ifelse(title2 == "" ,
                       str_trim(str_split(title,"\\]")[[1]][[3]]), title3 )) %>%
  ungroup() %>%
  select(-title, -title2) %>%
  rename(platform2=title3) %>%
  mutate(platform2=str_replace_all(platform2, " expression beadchip", "")) %>%
  mutate(platform2=tolower(str_replace_all(platform2, "_", " ")) ) %>%
  group_by(platform2) %>%
  mutate(platform2=str_split(platform2, " custom| based")[[1]][[1]])

gpl_info5 <- gpl_info4 %>%
  mutate(platform2=str_replace_all(platform2, "ver", "v")) %>%
  mutate(platform2=str_replace_all(platform2, "vsion ", "v")) %>%
  mutate(platform2=str_replace_all(platform2, "-", "")) %>%
  mutate(platform2=str_replace_all(platform2, "genechip", "")) %>%
  mutate(platform2=case_when(
    str_detect(platform2, "human genome u133 plus 2.0 array|hgu133 plus 2|133 plus 2.0|133 plus2|u133 v2") ~  
      "affymetrix human genome u133 plus 2.0 array",
    str_detect(platform2, "hgu133a" )~ 
                 "affymetrix human genome u133 array",
    str_detect(platform2, "hta2|hta 2| human transcriptome array 2.0") ~ 
      "affymetrix human genome hta 2.0 array",
    str_detect(platform2, "430 2.0 array") ~
      "affymetrix mouse genome 430 2.0 array",
    str_detect(platform2, "illumina humanwg6 v2") ~ "illumina humanwg6 v2.0" ,
    str_detect(platform2, "humanref8 v2") ~ 
                 "humanref8 v2.0",
    str_detect(platform2, "illumina mousewg6 v2") ~
      "illumina mousewg6 v2",
    str_detect(platform2, "illumina mouseref8 v2")  ~ "illumina mouseref8 v2",
    str_detect(platform2, "affymetrix human gene st 1.0") ~ "affymetrix human gene 1.0 st array",
    TRUE ~ platform2)) %>% 
  mutate(platform2=str_squish(str_trim(platform2))) %>%
  mutate(platform2=str_replace(platform2, " 0", "")) %>%
  mutate(platform2=str_replace(platform2, "ref| ref", "")) %>%
  mutate(platform2=str_replace(platform2, "\\.0", ""))

missing_src <- gpl_info5 %>% 
  filter(!str_detect(platform2, "illumina|affymetrix|sentrix|nanostring")) %>% 
  pull(gpl)
gpl_extra <- dbGetQuery(con, sprintf("SELECT gpl, title, manufacturer FROM gpl WHERE gpl IN ('%s');", 
                                     paste(missing_src, collapse="','")))
dbDisconnect(con)


missing_manufacturer <- gpl_info5 %>% 
  filter(!str_detect(platform2, "illumina|affymetrix|sentrix|nanostring")) %>%
  left_join(gpl_extra %>% select(gpl, manufacturer)) %>%
  mutate(manufacturer=tolower(manufacturer)) %>%
  mutate(manufacturer=case_when(
    str_detect(manufacturer, "illumina") ~ "illumina",
    str_detect(manufacturer, "affymetrix") ~ "affymetrix",
    TRUE ~ "")) %>%
  mutate(platform2=str_trim(paste(manufacturer, platform2, sep=" "))) %>%
  select(-manufacturer) %>% 
  ungroup()

gpl_info6 <- gpl_info5 %>% anti_join(missing_manufacturer, by="gpl") %>%
  bind_rows(missing_manufacturer) %>%
  select(-short_title) %>%
  mutate(platform2=ifelse(platform2=="human8 v2", "illumina human8 v2", platform2))


sample_metadata_w_plat <- sample_metadata2 %>% 
  left_join(gsm_to_gpl, by=c("sample_acc"="gsm")) %>%
  left_join(gpl_info6 %>% distinct(gpl, platform2))

other_plat <- sample_metadata_w_plat %>% 
  filter(is.na(platform2))  %>% 
  distinct(gpl, platform) %>% 
  group_by(platform) %>%
  mutate(platform3=tolower(str_trim(str_split(platform, "\\(")[[1]][[1]]))) %>%
  mutate(platform3=str_replace_all(platform3, "_", " "))  %>%
  mutate(platform3=str_replace_all(platform3, "version ", "v"))  %>%
  mutate(platform3=str_replace_all(platform3, "ref|ref ", "")) %>%
  mutate(platform3=str_replace_all(platform3,"-|\\.0", "")) %>%
  select(-gpl)

sample_metadata_w_plat2 <- sample_metadata_w_plat %>% 
  left_join(other_plat, by="platform") %>%
  mutate(cond_platform=ifelse(is.na(platform2), platform3, platform2))

sample_metadata_w_plat3 <- sample_metadata_w_plat2 %>% 
  select(-platform, -platform2, -platform3, -gpl) %>%
  mutate(platform=cond_platform) %>%
  select(-cond_platform) %>%
  distinct()

list_plat_final <- sample_metadata_w_plat3 %>% distinct(platform) %>% arrange(platform) 

list_plat_final %>% pull(platform)


#-------------- #
# deal with DUPLICATES
dup_samp <- sample_metadata_w_plat3 %>% 
  group_by(sample_acc) %>% 
  mutate(n=n()) %>% 
  filter(n>1) %>%
  arrange(sample_acc)
#dup_samp %>% select(-data_type) %>% distinct() %>% filter(duplicated(sample_acc))
# all of the duplicated are in BOTH microarray and RNA-seq

sample_metadata4 <- sample_metadata_w_plat3 %>% 
  semi_join(dup_samp %>% distinct(sample_acc)) %>%
              group_by(sample_acc) %>%
              summarize(data_type=paste(data_type, collapse=";"),
                        platform=unique(platform),
                        src=unique(src),
                        organism=unique(organism),
                        scrna=any(scrna)) %>%
  bind_rows(sample_metadata_w_plat3 %>% 
              anti_join(dup_samp %>% distinct(sample_acc)))
sample_metadata4 <- sample_metadata4 %>% ungroup() %>% mutate(sample_acc=str_trim(sample_acc))

study_metadata3 <- study_metadata2 %>%  
  group_by(study_acc) %>%
  summarize(data_type=paste(data_type, collapse=";"),
            title=paste(unique(title), collapse=";"),
            description=paste(unique(description), collapse=";"),
            date=paste(unique(date[date!="" & !is.na(date)]), collapse=";"),
            src=paste(unique(src),collapse=";"),
            organism=paste(unique(organism), collapse=";"),
            scrna=any(scrna))

# check assumptions!
stopifnot(length(unique(sample_metadata4$sample_acc))==nrow(sample_metadata4))
stopifnot((study_metadata3 %>% distinct(study_acc) %>% nrow())==nrow(study_metadata3))


# filter the rest
sample_to_study2 <- sample_to_study %>% 
  semi_join(study_metadata3, by="study_acc") %>%
  semi_join(sample_metadata4, by="sample_acc")
study_metadata4 <- study_metadata3 %>% semi_join(sample_to_study2, by="study_acc")
print(nrow(study_metadata3))
print(nrow(study_metadata4))
sample_metadata5 <- sample_metadata4 %>% semi_join(sample_to_study2, by="sample_acc")
print(nrow(sample_metadata4))
print(nrow(sample_metadata5))

stopifnot(length(setdiff(study_metadata4$study_acc, sample_to_study2$study_acc))==0)
stopifnot(length(setdiff(sample_metadata5$sample_acc, sample_to_study2$sample_acc))==0)

# write it out
list_plat_final %>% 
  write_csv("data/01_sample_lists/list_platforms.csv")

sample_metadata5 %>%
  write_csv("data/01_sample_lists/list_samples.csv")

study_metadata4 %>%
  write_csv("data/01_sample_lists/list_studies.csv")

sample_to_study2 %>%
  write_csv("data/01_sample_lists/sample_to_study.csv")

# TODO: remove mirna + tiling platforms?


# ---- patch missing metadata? ---- #

# TODO - grab date for 401 studies missing this info :/ 
# what ENA metadata do we need to download?
# -- run this

# what sample metadata are we missing?
