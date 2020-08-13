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

require('tidyverse')
require('data.table')
require('scales')  

#  --- set up color scheme --- #
library(RColorBrewer)
my.l <- brewer.pal(n = 8, name = 'Set2')
blues <- brewer.pal(9,name="Blues")
oranges <- brewer.pal(9,name="Oranges")
my.cols3 <- c (my.l[3],my.l[2], my.l[8])
my.cols4 <- c (my.l[3],my.l[4],my.l[2], my.l[8])
my.cols6 <- c (my.l[3], blues[4],my.l[4],oranges[4], my.l[2], my.l[8])



# ---- 2. CONSTRUCT DATASETS ---- #
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


# make this table 
#   | sample_acc | organism | data_type | platform | study_acc | metadata_sex | expr_sex

comb_metadata <- do.call(rbind, list(human_microarray, 
                                        human_rnaseq, 
                                        mouse_microarray, 
                                        mouse_rnaseq)) %>%
  select(sample_acc, organism, data_type, platform, study_acc, 
         metadata_sex, expr_sex, p_male ) 

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

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")
# ----- 3. COUNT TABLES ----- #
# TODO: UPDATE THE COUNTS  / overlap

# // TODO: what is our cutoff ?
cutoff <- comb_metadata %>%
  mutate(expr_sex=case_when(
    is.na(expr_sex) ~ "unknown",
    p_male < 0.7 & p_male > 0.3 ~ "unclear",
    TRUE ~ expr_sex
  ))


ggplot(cutoff, 
       aes(x=data_type, fill=expr_sex)) +
  geom_bar()+
  scale_fill_manual(values=my.cols4) +
  facet_grid(.~organism)


by_sample <- comb_metadata %>% 
  mutate(expr_sex=case_when(
    is.na(expr_sex) ~ "missing",
    p_male < 0.7 & p_male > 0.3 ~ "unlabeled",
    TRUE ~ expr_sex
  )) %>%
  mutate(metadata_sex=ifelse(is.na(metadata_sex) | metadata_sex=="unknown", 
                             "unlabeled", metadata_sex)) %>%
  filter(data_type=="microarray" |
           present & num_reads > 100000) %>%
  #select(-p_male) %>%
  pivot_longer(cols=c(expr_sex, metadata_sex), 
               names_to="label_type",
               values_to="sex_lab") %>%
  # change names so prettier
  mutate(label_type=case_when(label_type=="metadata_sex" ~ "metadata",
                              label_type=="expr_sex" ~ "expression")) %>%
  mutate(data_type=case_when(data_type=="RNA-seq" ~ "rnaseq",
                             TRUE ~ data_type)) %>%
  mutate(label_type=factor(label_type, levels=c("metadata", "expression")))


ggplot(by_sample %>%
         mutate(sex_lab = ifelse(sex_lab=="mixed", "mixed/pooled", sex_lab)) %>%
         mutate(sex_lab=factor(sex_lab, 
                               levels=c("female", "mixed/pooled", 
                                        "male", "unlabeled"))) %>%
         rename(`sample sex`=sex_lab), 
       aes(x=data_type, fill=`sample sex`))+
  geom_bar()+
  facet_grid(label_type ~organism)+  
  xlab("")+
  ylab("Number of samples")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(labels = comma)+
  scale_fill_manual(values=c(my.cols4) )

ggsave("figures/paper_figs/fig1_sample.png")

by_sample %>% write_csv("data/sample_metadata_filt.csv")

# by-study counts - NOTE: slow
by_study <- by_sample %>%
  separate_rows(study_acc, sep=";") %>%
  select(-present, -num_reads, -p_male) %>%
  group_by(study_acc, organism, data_type, label_type) %>%
  summarize(num_tot=n(),
         num_f=sum(sex_lab=="female", na.rm=TRUE),
         num_m=sum(sex_lab=="male", na.rm=TRUE),
         num_unlab=sum(sex_lab=="unlabeled", na.rm=TRUE))

by_study2 <- by_study %>%
  mutate(num_present=num_f+num_m) %>%
  mutate(study_sex=case_when(
    num_tot==num_f ~ "female only",
    num_tot==num_m ~ "male only",
    num_tot <= 60 & (0.5*num_tot > num_present) ~ "unknown",
    num_tot > 60 & (num_present < 30 ) ~ "unknown",
    num_f > 0.8*num_present ~ "mostly female",
    num_m > 0.8*num_present ~ "mostly male",
    TRUE ~ "mixed sex"
  )) %>%
  mutate(study_sex=factor(study_sex, levels=c("female only", "mostly female", 
                                              "mixed sex", "mostly male", 
                                              "male only", "unknown")))

by_study2 %>% write_csv("data/study_sex_lab.csv")
# SAVE THE BY-STUDY COUNTS!
#by_study2 <- read_csv("data/study_sex_lab.csv")
ggplot(by_study2 %>% 
         mutate(study_sex=factor(study_sex, levels=c("female only", "mostly female",
                                                     "mixed sex", "mostly male", 
                                                     "male only", "unknown"))) %>%
         mutate(label_type=factor(label_type, levels=c("metadata", "expression"))) %>%
         rename(`study sex`=study_sex), 
       aes(x=data_type, fill=`study sex`)) +
  geom_bar()+
  facet_grid(label_type~organism)+
  xlab("")+
  ylab("Number of studies")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=my.cols6)

  
ggsave("figures/paper_figs/fig1_study.png")

# --- create alluvial diagrams x 4 --- #
# microarray, rnaseq x human, mouse
require('ggalluvial')

flow_freq_counts <- by_sample %>%
    select(sample_acc, label_type, sex_lab, organism, data_type) %>%
    pivot_wider(id_cols=c(sample_acc,organism, data_type), names_from="label_type", values_from="sex_lab") %>%
    group_by(organism, data_type, metadata, expression) %>% 
    mutate(Freq=n()) %>% 
    select(-sample_acc) %>% 
    unique() %>%
    ungroup() %>%
    mutate(row_id=1:n()) %>%
    pivot_longer(cols=c("metadata", "expression"), 
                 names_to="labeling_method", values_to="sex") %>%
    mutate(sex=ifelse(sex %in% c("unknown","missing", "not labeled"), "unlabeled", sex)) %>%
    filter(sex != "mixed")   %>%
    mutate(row_id=as.factor(row_id), 
           labeling_method=factor(labeling_method, 
                                  levels=c("metadata", "expression")),
           sex=factor(sex, levels=c("female", "male", "unlabeled"))) %>%
    unique() 
  
  
ggplot(flow_freq_counts,
         aes(x = labeling_method, 
             stratum = sex, alluvium = row_id,
             y = Freq,
             fill = sex, label = sex)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3) +
    xlab("Label source")+ylab("Number of samples")+
    theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + 
    theme(legend.position = "none") +
    scale_fill_manual(values=my.cols3)+
    scale_y_continuous(labels = comma) +
    facet_grid(vars(data_type), vars(organism), scales="free")
ggsave("figures/paper_figs/fig1_sample_alluvial.png")  

# study level
study_flow_freq_counts <- by_study2 %>%
  mutate(study_sex=as.character(study_sex)) %>%
  mutate(study_sex=case_when( # for easier vis!
    (label_type=="metadata" & study_sex=="mostly male") ~ "mixed sex",
    (label_type=="metadata" & study_sex=="mostly female") ~ "mixed sex",
    TRUE ~ study_sex
    )) %>%
  select(study_acc, label_type, study_sex, organism, data_type) %>%
  pivot_wider(id_cols=c(study_acc,organism, data_type), names_from="label_type", 
              values_from="study_sex") %>%
  group_by(organism, data_type, metadata, expression) %>% 
  mutate(Freq=n()) %>% 
  select(-study_acc) %>% 
  unique() %>%
  ungroup() %>%
  mutate(row_id=1:n()) %>%
  pivot_longer(cols=c("metadata", "expression"), names_to="labeling_method", values_to="sex") %>%
  mutate(sex=ifelse(sex=="unknown", "unlabeled", sex)) %>%
  mutate(row_id=as.factor(row_id), 
         labeling_method=factor(labeling_method, levels=c("metadata", "expression")),
         sex=factor(sex, levels=c("female only", "mostly female",
                                     "mixed sex", "mostly male", 
                                     "male only", "unlabeled"))) %>%
  unique() 


ggplot(study_flow_freq_counts,
         aes(x = labeling_method, 
             stratum = sex, 
             alluvium = row_id,
             y = Freq,
             fill = sex, label = sex)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 2.5) +
    xlab("Label source")+ylab("Number of studies")+
    theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + 
    theme(legend.position = "none") +
    scale_fill_manual(values=my.cols6)+
    facet_grid(vars(data_type), vars(organism), scales="free")
  
ggsave("figures/paper_figs/fig1_study_alluvial.png")


# --- create Supplementary Table 1 --- #
# missingness - these tables
#  | organism | dataset | platform | sex | num_studies
#  | organism | dataset | platform | sex | num_samples
counts_sample <- by_sample %>%
  group_by(organism, data_type, label_type, sex_lab) %>%
  count() %>%
  ungroup() 

counts_sample_table <- counts_sample %>% 
  pivot_wider(names_from=sex_lab, values_from=n, 
              values_fill = list(n=0)) %>%
  mutate(num_samples=(female+male+mixed+unlabeled)) %>%
  mutate(frac_missing=unlabeled/num_samples,
         frac_female=female/num_samples,
         frac_male=male/num_samples) 


counts_sample_table %>%
  write_csv("tables/s1a_sample_sex_counts.csv")

counts_study <- by_study2  %>%
  group_by(organism, data_type, label_type, study_sex) %>%
  count() %>%
  ungroup()

counts_study_table <- counts_study %>% 
  group_by(organism, data_type, label_type) %>%
  mutate(num_studies=sum(n)) %>%
  ungroup() %>%
  pivot_wider(names_from=study_sex, values_from=n, values_fill = list(n=0)) %>%
  
  # // TODO: where does mostly go here?? should fraction be of labeled studies?
  mutate(frac_missing=unknown/num_studies,
         frac_male_only=`male only`/num_studies,
         frac_female_only=`female only`/num_studies,
         frac_mixed_sex=`mixed sex`/num_studies,
         frac_single_sex=(`male only`+`female only`)/num_studies) %>%
  select(-num_studies, -frac_missing, everything(), num_studies, frac_missing) 

# // TODO: format sig figs

counts_study_table %>%  
  write_csv("tables/s1b_study_sex_counts.csv")


# are the differences significant? can't do this! 
# b/c RNA-seq metadata is lossy
counts_sample_table %>% 
  filter(label_type=="metadata") %>% 
  select(organism, data_type, unlabeled, num_samples) %>%
  mutate(num_present=num_samples-unlabeled) %>%
  select(-num_samples) 

# # // TODO: stop hard-coding
# mat1 <- matrix(c(233526,96982,195258, 34531), byrow=TRUE, nrow=2)
# mat2 <- matrix(c(83926 ,39353, 298072,61070), byrow=TRUE, nrow=2)
# rownames(mat1) <- c("microarray", "rnaseq")
# colnames(mat1) <- c("missing", "present")
# human_sample_test <- chisq.test(mat1)
# 
# human_sample_test$expected
# human_sample_test$observed

# --- create Supplementary Table N --- #
# overlap


# <-- where did I look into this ("05_analysis/01_compendia_overlap.R")






