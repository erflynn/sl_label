# 02b_rnaseq_date_data.R
# E Flynn
# 07/20/2020
#
# Code for grabbing SRA date data
# This is complicated: 
#  - refine-bio does *NOT* have RNA-seq date data
#  - SRA has this in some places but it is unclear what field to use
# 
# One complicated piece of SRA data is that we are looking at studies ("SRP") 
# and runs ("SRR", but we call them samples), but there are submissions ("SRA"),
# samples ("SRS") and experiments ("SRX"). 
#
# There does appear to be a lot of cool information in the sample attribute data, but
# there is a lot of missingness. For now, we'll stick with the `updated_date`
# from SRA even though it is not ideal but we can be consistent.
#
# // TODO:
#  - check if BioProject would have this info and make it better?
#  - look at sample vs study
#  - use NER for dates?


# when we look in SRAdb, there are lots of missing data
# instead using ENA, only 18 missing

require('tidyverse')
require('rjson')
require('lubridate')

# the rb data has NO dates for RNA-seq
experiment_metadata <- read.csv("data/01_metadata/human_rnaseq_experiment_metadata.csv",
                                stringsAsFactors = FALSE) %>%
  bind_rows(read.csv("data/01_metadata/mouse_rnaseq_experiment_metadata.csv",
                     stringsAsFactors = FALSE) ) %>%
  select(study_acc, date) 
#THESE ARE ALL NAs

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")

rnaseq_studies <- (comb_metadata %>% filter(data_type=="rnaseq") %>%
                     select(study_acc) %>%
                     unique() %>%
                     separate_rows(study_acc, sep=";") %>%
                     unique())

rnaseq_studies %>% write_csv("data/rnaseq_studies.csv")


sra_dates <- fromJSON(file="data/dates/sra_date.json") 
sra_date_df <- do.call(rbind, lapply(sra_dates, function(x) data.frame(x, stringsAsFactors = FALSE)))
sra_date_df$accession <- rownames(sra_date_df)
table(sra_date_df$first_public=="") # none missing
sra_date_df %>% 
  select(accession, first_public) %>% 
  rename(date=first_public) %>%
  mutate(date=ymd(date)) %>%
  write_csv("data/dates/sra_study_ena.csv")

nrow(rnaseq_studies)-nrow(sra_date_df)
# "sra_missing_dates.txt" contains the 18 studies that are missing

# ---- wound up parsing server-side from ENA --- #

# ---- STOP HERE ---- #
require('SRAdb')



# -- use SRA db to grab data
sqlfile <- file.path("SRAmetadb.sqlite")
sra_con <- dbConnect(SQLite(),sqlfile)

date_data <- dbGetQuery(sra_con, 
                        "SELECT run_accession, run_date, study_accession, updated_date, submission_accession, submission_date FROM sra;")

save(date_data, file="sra_date.RData")
dbDisconnect(sra_con)

load("data/dates/sra_date.RData") # --> date_data

# -- problem: the `submission_date` and `run_date` fields have high missingness
list_rnaseq_samples <- (comb_metadata %>% filter(data_type=="rnaseq"))$sample_acc
list_rnaseq_studies <- (comb_metadata %>% filter(data_type=="rnaseq") %>% 
                          select(study_acc) %>% 
                          unique() %>%
                          separate_rows(study_acc, sep=";") %>%
                          unique())$study_acc

sra_date_samples <- date_data %>% 
  filter(run_accession %in% list_rnaseq_samples)
sra_date_studies <- date_data %>% 
  filter(study_accession %in% list_rnaseq_studies)

length(setdiff(list_rnaseq_samples, date_data$run_accession)) # 5k are missing :(
table(is.na(sra_date_samples$submission_date)) # most are missing dates


sra_date_studies2 <- sra_date_studies %>% 
  filter(!is.na(submission_date)) %>% 
  select(study_accession, updated_date, submission_accession, submission_date) %>% 
  unique()

study_sra <- sra_date_studies %>% 
  select(study_accession, updated_date, submission_date) %>%
  unique()
table(is.na(study_sra$updated_date))
table(is.na(study_sra$submission_date)) # mostly missing

study_sra %>% 
  select(study_accession, updated_date) %>%
  write_csv("data/dates/sra_study_date.csv")


# -- alternately - we could use `updated_date`, which is present in everything
# let's look at the relationship of updated_date to submission and run dates
date_comp_df <- sra_date_samples %>% 
  select(run_accession, run_date, updated_date, submission_date) %>%
  as_tibble() %>%
  mutate(run_date=ymd(run_date), 
         updated_date=ymd(updated_date),
         submission_date=ymd(submission_date)) %>%
  mutate(run_to_update=updated_date-run_date,
         run_to_submission=submission_date-run_date,
         submission_to_update=updated_date-submission_date)
# 9 failed to parse... hmm...
date_comp_df %>% filter(is.na(updated_date)) # all present
date_comp_df %>% 
  select(run_accession, updated_date) %>% 
  write_csv("data/dates/sra_date_info.csv")

ggplot(date_comp_df, aes(x=run_to_update/365))+geom_histogram() 
# ^ some of these have big gaps, most within 5y

ggplot(date_comp_df, aes(x=submission_to_update/365))+geom_histogram()
# more of these are w/in 5y, tho mode is ~3

ggplot(date_comp_df, aes(x=run_to_submission/365))+geom_histogram()
# seems to be around 0.5y (mode), but variable

# it is confusing that the run-->update period be negative
date_comp_df %>% filter(run_to_update < 0) # dates in the future?

# so: update present, but mb not what we want?

# write this out!

# -- let's try looking at the attribute data for samples, etc --- #
sra_con <- dbConnect(SQLite(), sqlfile)
h_exp_sample <- read_csv("refine_bio/data/01_metadata/human_rnaseq_exp_to_sample.csv")
m_exp_sample <- read_csv("refine_bio/data/01_metadata/mouse_rnaseq_exp_to_sample.csv")
h_exp_sample <- read_csv("data/01_metadata/human_rnaseq_exp_to_sample.csv")
m_exp_sample <- read_csv("data/01_metadata/mouse_rnaseq_exp_to_sample.csv")



# convert runs to SRA samples
h_runs_to_sample <- sraConvert(h_exp_sample$sample_acc, 'sample', sra_con) 
m_runs_to_sample <- sraConvert(m_exp_sample$sample_acc, 'sample', sra_con) 
save(h_runs_to_sample, m_runs_to_sample, file="run_sample_sra.RData")

# sample attributes -- abt half have something (tbd if dates tho!)
sample_list <- unique(c(h_runs_to_sample$sample, m_runs_to_sample$sample))
sample_list_str <- paste(sample_list, collapse="\',\'")
sample_ds <- dbGetQuery(sra_con, 
                        sprintf("SELECT sample_accession, sample_attribute FROM sample WHERE sample_accession IN ('%s')", 
                                sample_list_str))
save(sample_ds, file="sample_attr_sra.RData")
table(is.na(sample_ds$sample_attribute))
#FALSE   TRUE 
# 273638 250724 


# run attributes -- tho mostly empty
run_list <- unique(c(h_runs_to_sample$run, m_runs_to_sample$run))
run_list_str <- paste(run_list, collapse="\',\'")
run_ds <- dbGetQuery(sra_con, 
                     sprintf("SELECT run_accession, run_attribute FROM run WHERE run_accession IN ('%s')", 
                             run_list_str))
save(run_ds, file="run_attr_sra.RData")

table(is.na(run_ds$run_attribute))
# FALSE   TRUE 
# 133062 568748 

# study attributes (only 1/3 have something)
study_list <- unique(c(h_exp_sample$study_acc, m_exp_sample$study_acc))
study_list_str <- paste(study_list, collapse="\',\'")
study_ds <- dbGetQuery(sra_con, sprintf("SELECT study_accession, study_attribute FROM study WHERE study_accession IN ('%s')", 
                                        study_list_str))
save(study_ds, file="study_attr_sra.RData")


table(is.na(study_ds$study_attribute))
# FALSE  TRUE 
# 3438  9239

# skipping experiment?
#dbGetQuery(sra_con,"SELECT experiment_accession, experiment_attribute FROM experiment WHERE experiment_attribute IS NOT NULL LIMIT 50;")

dbDisconnect(sra_con)

# -- try parsing the sample attributes for dates -- #
load("data/dates/sample_attr_sra.RData") # --> sample_ds

ds <- sample_ds %>%
  filter(!is.na(sample_attribute)) %>%
  as_tibble() %>%
  #sample_n(300) %>% # did this for randomization, TODO remove and run on everything
  separate_rows(sample_attribute, sep=" \\|\\| ") %>%
  separate(sample_attribute, into=c("key", "value"), sep=": ", extra="merge")

unique(ds$key)
ds %>% group_by(key) %>% count() %>% arrange(desc(n))
keys <- (ds %>% group_by(key) %>% count() %>% arrange(desc(n)))$key

keys[str_detect(sapply(keys,tolower), "sex|gender")]
date_keys <- keys[str_detect(sapply(keys,tolower), "date")] 
# ones to exclude:
filt_date_keys <- date_keys[!str_detect(sapply(date_keys,tolower), "birth|death|biopsy|diagnosis")]
subm_keys <- keys[str_detect(sapply(keys,tolower), "subm|pub")] 
# keep: "INDSC first public", "ENA-FIRST-PUBLIC"

length(unique(sample_ds$sample_accession)) # 524362
ds %>% filter(key %in% c(filt_date_keys,"INSDC first public", "ENA-FIRST-PUBLIC"))  # 8347
#%>% select(sample_accession) %>% unique() # 5011

ds_date <- ds %>% filter(key %in% c(filt_date_keys,"INSDC first public", "ENA-FIRST-PUBLIC")) %>% 
  mutate(val_date=ymd(value), val_date2=dmy(value), val_date3=mdy(value),
         val_date4=ymd_hms(value)) 
ds_date %>%
  filter(is.na(val_date) & is.na(val_date2) & is.na(val_date3) & is.na(val_date4)) %>%
  select(value) %>% unique() 
# this grabs pretty much all of them -- yay!
# BUT we need to decide a) which is the consensus date and b) whether we want it

# check if there are fields we are missing
#extra_dates <- ds %>% filter(!key %in% c(filt_date_keys,"INSDC first public", "ENA-FIRST-PUBLIC")) %>%
#  mutate(val_date=ymd(value), val_date2=dmy(value), val_date3=mdy(value),
#             val_date4=ymd_hms(value)) %>%
#  filter(!(is.na(val_date) & is.na(val_date2) & is.na(val_date3) & is.na(val_date4)))
#extra_dates %>% select(-sample_accession) %>% unique() %>% View()
# THESE ARE NOT DATES. skip

# // TODO: we really need good NER for these data!

# Summarizing what we did:
#  there *IS* some date data in the samples, etc

load("data/dates/study_attr_sra.RData") # study_ds

# 547
study_ds %>% filter(!is.na(study_attribute)) %>% filter(str_detect(study_attribute, "ENA-FIRST-PUBLIC")) %>% nrow()
table(is.na(study_ds$study_attribute))
study_ds %>% filter(is.na(study_attribute)) %>% sample_n(10)

# if we look online -- this is MUCH better