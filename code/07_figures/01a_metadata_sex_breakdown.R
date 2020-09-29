# 01a_metadata_sex_breakdown.R
# E Flynn
# 7/20/2020
# Updated 8/12/2020 to add in missing RNA-seq samples
# 
# Code for looking at the metadata sex breakdown. 
# This generates by-study and by-sample tables used for follow up analyses.
# We also plot the alluvial diagrams.
#
# TODO: 
# - when to exclude problem platforms?
# - separate out into multiple files
# CONCERN: MOUSE RNA-seq shows a high fraction of female samples labeled ambiguous
#  - this indicates model misspecification!

require('tidyverse')
require('data.table')
require('scales')  



# ---- 1. CONSTRUCT DATASETS ---- #
constructDS <- function(organism, data_type, filter_samples){
  # metadata sex labels
  metadata_sex <- fread(sprintf("data/01_metadata/%s_%s_metadata_sex.csv", 
                                organism, data_type),
        data.table=FALSE, stringsAsFactors=FALSE) %>%
    as_tibble() %>%
    select(-sex) %>%
    rename(sample_acc=acc, metadata_sex=mapped_sex) %>%
    unique()

  stopifnot(length(unique(metadata_sex$sample_acc))==nrow(metadata_sex))
  
  # filter to remove specific samples -- we do this for microarray, remove rnaseq
  if (!is.null(filter_samples) & data_type=="microarray"){
    metadata_sex <- metadata_sex %>% 
      filter(!sample_acc %in% filter_samples)
  }
  print("Metadata sex")
  
  # add studies
  sample_str <- ifelse(data_type=="microarray", "", "_rnaseq")
  exp_to_sample <- read_csv(sprintf("data/01_metadata/%s%s_exp_to_sample.csv", organism, sample_str))
  sample_to_study <- exp_to_sample %>% 
    group_by(sample_acc) %>% 
    summarize(study_acc=paste(unique(study_acc), collapse=";")) %>%
    ungroup()
  metadata_sex_w_study <- metadata_sex %>%
    left_join(sample_to_study, by=c("sample_acc"))
  stopifnot(nrow(metadata_sex_w_study)==nrow(metadata_sex))
  
  print("Studies")
  
  
  # add platform + imputed sex
  if (data_type=="microarray"){
    plat_dat <- fread(sprintf("data/01_metadata/%s_metadata.csv", organism), # platform
                data.table=FALSE, stringsAsFactors=FALSE) %>% 
      as_tibble() %>%
      rename(sample_acc=acc)
    sl_dat <- fread(sprintf("data/%s_metadata2.csv", organism), # sex_lab, pred
                   data.table=FALSE, stringsAsFactors=FALSE) %>% 
      as_tibble() %>%
      rename(sample_acc=acc)
  } else {
    plat_dat <- fread(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", organism),
                data.table=FALSE, stringsAsFactors=FALSE) %>% 
      as_tibble() %>%
      rename(sample_acc=acc)
    sl_dat <- fread(sprintf("data/02_labeled_data/%s_rnaseq_sl.csv", organism),
                   data.table=FALSE, stringsAsFactors=FALSE) %>% 
      as_tibble() %>%
      rename(sample_acc=id) %>%
      mutate(sex_lab=ifelse(pred > 0.5, "male", "female")) 
    
  }
  metadata_w_plat <- metadata_sex_w_study %>% 
    left_join(plat_dat %>% 
                select(sample_acc, platform) %>%
                unique())
  stopifnot(nrow(metadata_sex_w_study)==nrow(metadata_w_plat))
  print("Platform")
  metadata_w_sl <- metadata_w_plat %>% 
    left_join(sl_dat %>% 
                select(sample_acc, sex_lab, pred) %>% 
                rename(expr_sex=sex_lab, p_male=pred) %>%
                unique())
  stopifnot(nrow(metadata_sex_w_study)==nrow(metadata_w_sl))
  
  # add labels
  metadata_w_sl$organism <- organism
  metadata_w_sl$data_type <- data_type
  return(metadata_w_sl)
}

mouse_rnaseq <- constructDS("mouse", "rnaseq")
stopifnot(length(unique(mouse_rnaseq$sample_acc))==nrow(mouse_rnaseq))
mouse_microarray <- constructDS("mouse", "microarray", mouse_rnaseq$sample_acc)
human_rnaseq <- constructDS("human", "rnaseq")
human_microarray <- constructDS("human", "microarray", human_rnaseq$sample_acc)

# re-read all the sex labeling metadata
rnaseq_files <- list.files(path="data/08_model_out/", pattern="*.csv")
rnaseq_sl <- rnaseq_files %>% map(~read_csv(sprintf("data/08_model_out/%s",.x))) %>% bind_rows()
head(rnaseq_sl)




# make this table 
#   | sample_acc | organism | data_type | platform | study_acc | metadata_sex | expr_sex

comb_metadata <- do.call(rbind, list(human_microarray, 
                                        human_rnaseq, 
                                        mouse_microarray, 
                                        mouse_rnaseq)) %>%
  select(sample_acc, organism, data_type, platform, study_acc, 
         metadata_sex, expr_sex, p_male ) 

rnaseq_sl2 <- comb_metadata %>% filter(data_type=="rnaseq") %>% 
  full_join(rnaseq_sl,by=c("sample_acc"="id"))
ggplot(rnaseq_sl2 %>% filter(num_reads > 100000), aes(x=pred, col=metadata_sex))+
  geom_density()+
  facet_grid(organism ~ ., scales="free")

df <- rnaseq_sl2 %>% 
  filter(metadata_sex %in% c("male", "female"), 
         num_reads > 100000,
        !is.na(pred)) %>%
  mutate(est_sex=ifelse(pred > 0.5, "male", "female")) %>%
  filter(organism=="human")
table(df$est_sex==df$metadata_sex)/nrow(df)
table(df$expr_sex==df$metadata_sex)/nrow(df)
ggplot(df %>% filter(num_reads > 100000), aes(x=pred,  col=metadata_sex))+
  geom_density()+
  facet_grid(organism ~ ., scales="free")

# ---- read in QC information for RNA-seq ---- #
h_present <- read_csv("data/01_metadata/human_rnaseq_exp_to_sample2.csv") %>% mutate(organism="human")
m_present <- read_csv("data/01_metadata/mouse_rnaseq_exp_to_sample2.csv") %>% mutate(organism="mouse")

h_es <- read_csv("data/01_metadata/human_exp_to_sample_counts.csv") %>% select(-present, -study_acc)
m_es <- read_csv("data/01_metadata/mouse_exp_to_sample_counts.csv") %>% select(-present, -study_acc)
qc_data <- h_present %>% left_join(h_es) %>% 
  bind_rows(m_present %>% left_join(m_es) ) %>%
  select(organism, study_acc, sample_acc, present, num_reads)

missing_read_cts <- qc_data %>% filter(present & is.na(num_reads)) 
missing_read_cts %>% write_csv("data/missing_read_counts.csv")

# missing read cts
h_m <- read_csv("data/human_read_counts_m.csv") %>% 
  mutate(organism="human") %>% 
  select(organism, everything())
m_m <- read_csv("data/mouse_read_counts_m.csv") %>%
  mutate(organism="mouse") %>% 
  select(organism, everything())
fill_in_rc <- h_m %>% bind_rows(m_m)

# this fills in it all!
missing_read_cts %>% select(-num_reads) %>%
  left_join(fill_in_rc %>% 
              select(sample_acc, num_reads)) 

qc_data2 <- qc_data %>% 
  anti_join(fill_in_rc, by="sample_acc") %>%
  bind_rows(fill_in_rc)

comb_metadata_w_qc <- comb_metadata %>% 
  left_join(qc_data2)

comb_metadata_w_qc %>% 
  filter(!is.na(num_reads) & num_reads < 100000)  # removes 6.8k


comb_metadata_w_qc %>%
  filter(data_type=="rnaseq" & is.na(p_male) & present &
          num_reads > 100000) %>%
  select(-platform, -metadata_sex, -expr_sex) %>%
  select(organism, study_acc, sample_acc, present, num_reads) %>%
  write_csv("data/missing_rnaseq_pred.csv")
  #%>% # 4k
  #anti_join(missing_read_cts, by="sample_acc") # 588 

missing_dat <- read_csv("data/human_rnaseq_sl_m.csv") %>%
  bind_rows( read_csv("data/mouse_rnaseq_sl_m.csv") )

missing_filled_in <- comb_metadata_w_qc %>%
  filter(data_type=="rnaseq" & is.na(p_male) & present &
           num_reads > 100000) %>%
  left_join(missing_dat, by=c("sample_acc"="id")) %>%
  select(-p_male, -expr_sex) %>%
  rename(p_male=pred) %>%
  mutate(expr_sex=ifelse(p_male > 0.5, "male", "female")) %>%
  select(colnames(comb_metadata_w_qc))

comb_metadata_w_qc2 <- comb_metadata_w_qc %>%
  anti_join(missing_filled_in, by="sample_acc") %>%
  bind_rows(missing_filled_in)

stopifnot(length(unique(comb_metadata_w_qc2$sample_acc))==nrow(comb_metadata_w_qc2))
write_csv(comb_metadata_w_qc2, "data/01_metadata/combined_human_mouse_meta.csv")


comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", col_types="cccccccdld")
sl_updated <- read_csv("data/01_metadata/mapped_sl_all.csv")

table(sl_updated$mapped_sex)
table(comb_metadata$metadata_sex)
comb_metadata2 <- comb_metadata %>% 
  left_join(sl_updated %>% distinct(sample_acc, mapped_sex)) %>%
  mutate(mapped_sex=ifelse(is.na(mapped_sex), "unknown", mapped_sex))
stopifnot(nrow(comb_metadata)==nrow(comb_metadata2))

table(comb_metadata2$metadata_sex, comb_metadata2$mapped_sex)
mismatch <- comb_metadata2 %>% filter(mapped_sex %in% c("male", "female", "mixed") & 
                            metadata_sex %in% c("male", "female", "mixed") & 
                            metadata_sex != mapped_sex) 
mismatch %>% View()
#"GSE55231;GSE55232"
#"GSE73461;GSE73464"
#"GSE83951"
#"GSE1919"
#"GSE59512"
#"GSE65038;GSE65048"
sl_updated %>% filter(sample_acc=="GSM1332058")


comb_metadata2 %>% 
  group_by(data_type) %>%
  filter(mapped_sex!="unknown" & metadata_sex=="unknown") %>% count()
# 48353 (4.64%), mostly rescues RNA-seq


comb_metadata2 %>% 
  group_by(data_type) %>%
  filter(metadata_sex!="unknown" & mapped_sex=="unknown") %>% count() # 2302

comb_metadata3 <- comb_metadata2 %>% 
  select(-metadata_sex) %>% 
  rename(metadata_sex=mapped_sex)

comb_metadata3 %>% select(setdiff(colnames(comb_metadata2), "mapped_sex")) %>%
                            write_csv("data/01_metadata/combined_human_mouse_meta_v2.csv")

comb_metadata %>% filter(data_type=="rnaseq") %>% 
  distinct(sample_acc) %>% 
  write_csv("data/rnaseq_runs.csv", col_names=FALSE)



