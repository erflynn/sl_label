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


require('tidyverse')
require('SRAdb')
require('lubridate')

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")

# the rb data has NO dates for RNA-seq
experiment_metadata <- read.csv("data/01_metadata/human_rnaseq_experiment_metadata.csv",
                                stringsAsFactors = FALSE) %>%
  bind_rows(read.csv("data/01_metadata/mouse_rnaseq_experiment_metadata.csv",
                     stringsAsFactors = FALSE) ) %>%
  select(study_acc, date) 
#THESE ARE ALL NAs

# -- use SRA db to grab data
sqlfile <- file.path("SRAmetadb.sqlite")
sra_con <- dbConnect(SQLite(),sqlfile)

date_data <- dbGetQuery(sra_con, 
                        "SELECT run_accession, run_date, study_accession, updated_date, submission_accession, submission_date FROM sra;")

save(date_data, file="sra_date.RData")
dbDisconnect(sra_con)

load("data/sra_date.RData") # --> date_data

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
table(is.na(sra_date_samples$submission_date)) # SHIT most are missing dates


sra_date_studies2 <- sra_date_studies %>% 
  filter(!is.na(submission_date)) %>% 
  select(study_accession, updated_date, submission_accession, submission_date) %>% 
  unique()


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

ggplot(date_comp_df, aes(x=run_to_update/365))+geom_histogram() 
# ^ some of these have big gaps, most within 5y

ggplot(date_comp_df, aes(x=submission_to_update/365))+geom_histogram()
# more of these are w/in 5y, tho mode is ~3

ggplot(date_comp_df, aes(x=run_to_submission/365))+geom_histogram()
# seems to be around 0.5y (mode), but variable

# it is confusing that the run-->update period be negative
date_comp_df %>% filter(run_to_update < 0) # dates in the future?

# so: update present, but mb not what we want?


# -- let's try looking at the attribute data for samples, etc --- #
sra_con <- dbConnect(SQLite(), sqlfile)
h_exp_sample <- read_csv("refine_bio/data/01_metadata/human_rnaseq_exp_to_sample.csv")
m_exp_sample <- read_csv("refine_bio/data/01_metadata/mouse_rnaseq_exp_to_sample.csv")

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
ds <- sample_ds %>%
  filter(!is.na(sample_attribute)) %>%
  as_tibble() %>%
  sample_n(300) %>% # did this for randomization, TODO remove and run on everything
  separate_rows(sample_attribute, sep=" \\|\\| ") %>%
  separate(sample_attribute, into=c("key", "value"), sep=": ", extra="merge")

unique(ds$key)
ds %>% group_by(key) %>% count() %>% arrange(desc(n))
keys <- (ds %>% group_by(key) %>% count() %>% arrange(desc(n)))$key

keys[str_detect(sapply(keys,tolower), "sex|gender")]
date_keys <- keys[str_detect(keys, "date")] # this gives us 5 dates... (Adding time doesnt help)
subm_keys <- keys[str_detect(keys, "subm|pub")] # gives us 2 more!


