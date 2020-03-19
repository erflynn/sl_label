require('tidyverse')
require('lubridate')
options(stringsAsFactors = FALSE)

# read in metadata

get_dates <- function(prefix){
  experiment_metadata <- read.csv(sprintf("data/01_metadata/%s_experiment_metadata.csv",prefix))
  #table(is.na(experiment_metadata$date))
  #table(experiment_metadata$date!="") # we only have abt half
  
  experiment_metadata %>% 
    filter(date!="") %>%
    group_by(study_acc) %>%
    mutate(sm_date=strsplit(date, "T")[[1]][[1]])  %>%
    mutate(formatted_date=ymd(sm_date)) %>%
    select(study_acc, formatted_date) %>%
    as_tibble() %>% 
    ungroup()
}

reform_summary <- function(summ){
  summ %>%
    filter(labeling_method =="exprsex") %>% 
    select(study, sex) %>%
    mutate(sex=as.character(sex)) %>%
    mutate(study_sex=ifelse(sex %in% c("mostly-male", "mostly-female"), "mixed", sex)) %>%
    select(-sex) %>%
    rename(study_acc=study)
}

reform_summary_metadata <- function(summ){
  summ %>%
    filter(labeling_method =="metadata") %>% 
    select(study, sex) %>%
    mutate(sex=as.character(sex)) %>%
    mutate(study_sex=ifelse(sex %in% c("mostly-male", "mostly-female"), "mixed", sex)) %>%
    select(-sex) %>%
    rename(study_acc=study)
}


# read in study info
load("data/summary_files.RData") # --> h_summ
load("data/sample_level_sex.RData")
h_summ2 <- reform_summary(h_summ)
h_summ2_m <- reform_summary_metadata(h_summ)
h_dat <- get_dates("human") 
h_dat_plus <- inner_join(h_dat, h_summ2) %>%
  mutate(organism="human")
h_dat_plus_m <- inner_join(h_dat, h_summ2_m) %>%
  mutate(organism="human")


m_summ2 <- reform_summary(m_summ)
m_summ2_m <- reform_summary_metadata(m_summ)
m_dat <- get_dates("mouse") 
m_dat_plus <- inner_join(m_dat, m_summ2) %>%
  mutate(organism="mouse")
m_dat_plus_m <- inner_join(m_dat, m_summ2_m) %>%
  mutate(organism="mouse")

r_summ2 <- reform_summary(r_summ)
r_summ2_m <- reform_summary_metadata(r_summ)
r_dat <- get_dates("rat") 
r_dat_plus <- inner_join(r_dat, r_summ2) %>%
  mutate(organism="rat")
r_dat_plus_m <- inner_join(r_dat, r_summ2_m) %>%
  mutate(organism="rat")

dat_plus <- do.call(rbind, list(h_dat_plus, m_dat_plus, r_dat_plus)) %>%
  mutate(study_sex=factor(study_sex))

dat_plus_m <- do.call(rbind, list(h_dat_plus_m, m_dat_plus_m, r_dat_plus_m)) %>%
  mutate(study_sex=factor(study_sex))

ggplot(dat_plus, aes(x=formatted_date))+
  geom_histogram(aes( fill=study_sex))+
  ylab("number of studies")+
  xlab("date submitted")+
  facet_wrap(.  ~organism, ncol=1, scales="free")

ggplot(dat_plus_m, aes(x=formatted_date))+
  geom_histogram(aes( fill=study_sex))+
  ylab("number of studies")+
  xlab("date submitted")+
  facet_wrap(.  ~organism, ncol=1, scales="free")

dat_plus_both <- rbind(
  dat_plus_m %>% mutate(source="metadata"), 
  dat_plus %>% mutate(source="expression")
) %>%
  mutate(source=factor(source, levels=c("metadata", "expression")),
         study_sex=factor(study_sex, levels=c("female-only", "mixed", "male-only", "unlabeled")))

ggplot(dat_plus_both, aes(x=formatted_date))+
  geom_histogram(aes( fill=study_sex), bins=20)+
  ylab("number of studies")+
  xlab("date submitted")+
  facet_grid(rows=vars(organism), cols=vars(source),  scales="free")+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())  


# fraction single sex
require('zoo')
dat_plus2 <- dat_plus %>%
  arrange(organism, study_sex,formatted_date) %>%
  mutate(value=1) %>%
  group_by(organism,study_sex) %>%
  rename(date=formatted_date) %>%
  complete(date = full_seq(date, period = 1), fill = list(value = 0)) %>%
  group_by(organism, study_sex, month3=floor_date(date, "3 months")) %>%
  summarize(value=sum(value)) %>%
  mutate(tot= rollapplyr(value, width = 4, FUN = sum, partial = TRUE)) %>%
  select(organism, study_sex, month3, tot) %>%
  ungroup() %>%
  unique() %>%
  pivot_wider(id_cols=c(organism, month3), names_from=study_sex, values_from=tot, values_fn=list(tot=max)) %>%
  group_by(organism, month3) %>%
  mutate(frac_single_sex=(`female-only`+`male-only`)/(mixed+1),
         frac_m_only=`male-only`/(`female-only`+`male-only`+1))

# fraction male only and single sex
ggplot(dat_plus2, aes(x=month3, y=frac_single_sex))+geom_line(aes(col=organism))
ggplot(dat_plus2, aes(x=month3, y=frac_m_only))+geom_line(aes(col=organism))
  
# fraction unlabeled
#human_sl <- read_csv("data/02_labeled_data/human_all_sl.csv")

metadata <- read.csv("data/01_metadata/human_metadata.csv") # BE CAREFUL
metadata_filt <- metadata %>% left_join(exp_to_sample, by=c("acc"="sample_acc")) %>%
  semi_join(h_summ2, by=c("study_acc"))
no_cl <- metadata %>% filter(cl_line=="" )
exp_to_sample <- read.csv("data/01_metadata/human_exp_to_sample.csv")
exp2 <- exp_to_sample %>% filter(sample_acc %in% no_cl$acc)
h_dat_plus %>% 
  mutate(cl=ifelse(study_acc %in% exp2$study_acc, "tissue", "cell line")) %>%
ggplot( aes(x=formatted_date))+
    geom_histogram(aes( fill=study_sex))+
    ylab("number of studies")+
    xlab("date submitted") +
  facet_wrap(.~cl, ncol=1, scales="free")+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())  

# other plots:
# [x] how does this vary by organism
# [x] does the fraction unlabeled vary?
# [x] does this vary by cell vs tissue data?
#  -- are there cell line trends?
# 4. what's the deal w RNA-seq data?

