
# END GOAL:
# have these tables so we can make a figure:
# 1. sample | organism | data_type | sex | studies | sample_date
# 2. study | study_date
require('tidyverse')

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")
# --- we need to add date data... --- #

# 1. get GEOmetadb
require('GEOmetadb')
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

con <- dbConnect(SQLite(), "../labeling/GEOmetadb.sqlite") # sept 23, 2019
list_studies_str <- paste(list_studies, collapse="\',\'")
dbListTables(con)
res <- dbGetQuery(con, sprintf("SELECT gse.gse, submission_date FROM gse WHERE gse.gse IN ('%s');", list_studies_str))

stopifnot(nrow(res)==length(unique(res$gse)))
stopifnot(length(list_studies[str_detect(list_studies, "GSE")])==nrow(res))

# no missing dates! woot
any(is.na(res$submission_date))
any(res$submission_date=="")

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
table(str_detect(list_samples, "GSM")) #5k samples
length(list_samples[str_detect(list_samples, "GSM")])

dbDisconnect(con)
# missing 10 GSM samples


# 3. Grab ArrayExpress metadata
# https://www.ebi.ac.uk/arrayexpress/help/programmatic_access.html#Experiment_details
# > https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/E-xxxx-nnnnn
# > https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/E-xxxx-nnnnn/samples
# this is XML, so will have to download and parse


