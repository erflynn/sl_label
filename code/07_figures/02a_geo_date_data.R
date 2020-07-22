# 02a_geo_date_data.R
# E Flynn
# 07/20/2020
#
# Get data for microarray samples. Some of this is in refine-bio studies but
# it has high missingness --> we will query GEOmetadb and ArrayExpress directly.
#
# END GOAL:
# have these tables so we can make a figure:
# 1. sample | organism | data_type | sex | studies | sample_date
# 2. study | study_date

# -- UPDATE: we don't need to query GEOmetadb for this info! ignore. -- #

require('tidyverse')
require('GEOmetadb')


comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")

# 1. get GEO studies
# study-level
list_studies <- (comb_metadata %>% 
                   filter(data_type=="microarray") %>%
  select(study_acc) %>% 
  unique() %>%
  separate_rows(study_acc, sep=";") %>%
  unique())$study_acc

table(str_detect(list_studies, "GSE"))
# FALSE  TRUE 
# 260  19963 
non_gse_studies <- list_studies[!str_detect(list_studies, "GSE")]
non_gse_studies %>% as_tibble() %>% rename(study=value) %>% write_csv("data/dates/non_gse.csv")


con <- dbConnect(SQLite(), "../labeling/GEOmetadb.sqlite") # sept 23, 2019
list_studies_str <- paste(list_studies, collapse="\',\'")
dbListTables(con)
study_res <- dbGetQuery(con, sprintf("SELECT gse.gse, submission_date FROM gse WHERE gse.gse IN ('%s');", list_studies_str))

stopifnot(nrow(study_res)==length(unique(study_res$gse)))
stopifnot(length(list_studies[str_detect(list_studies, "GSE")])==nrow(study_res))

# no missing dates! woot
any(is.na(study_res$submission_date))
any(study_res$submission_date=="")

# sample-level
dbListFields(con, "gsm")
list_samples <- (comb_metadata %>% 
                   filter(data_type=="microarray"))$sample_acc
list_samples_str <- paste(list_samples, collapse="\',\'")

sample_res <- dbGetQuery(con, sprintf("SELECT gsm.gsm, submission_date FROM gsm WHERE gsm.gsm IN ('%s');", 
                                      list_samples_str))

any(is.na(sample_res$submission_date))
any(sample_res$submission_date=="")

length(list_samples)
nrow(sample_res)
table(str_detect(list_samples, "GSM")) #5k samples are non-GSM
length(list_samples[str_detect(list_samples, "GSM")])

dbDisconnect(con)
# missing 10 GSM samples from the list of dates, not bad

# SAVE THIS!
sample_res %>% write_csv("data/dates/geo_sample_dates.csv")
study_res %>% write_csv("data/dates/geo_study_dates.csv")

