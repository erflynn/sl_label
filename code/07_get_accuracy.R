
require('tidyverse')

human <- read.csv("data/human_metadata2.csv")
mouse <- read.csv("data/mouse_metadata2.csv")
rat <- read.csv("data/rat_metadata2.csv")


## get the study level labels 
h_map <- read.csv("data/human_exp_to_sample.csv")
m_map <- read.csv("data/mouse_exp_to_sample.csv")
r_map <- read.csv("data/rat_exp_to_sample.csv")

## clean up the metadata sex labels
# -- this is done manually into male, female, mixed, and unknown
# human %>% filter(!is.na(sex)) %>% select(sex) %>% unique() %>% write_csv("human_alt_sex_annot.csv")
# mouse %>% filter(!is.na(sex)) %>% select(sex) %>% unique() %>% write_csv("mouse_alt_sex_annot.csv")
# rat %>% filter(!is.na(sex)) %>% select(sex) %>% arrange(sex) %>% unique() %>% write_csv("rat_alt_sex_annot.csv")
mouse_alt <- read_csv("mouse_alt_sex_annot.csv")
human_alt <- read_csv("human_alt_sex_annot.csv")
rat_alt <- read_csv("rat_alt_sex_annot.csv")

# ------ reformat ----- #
addAltReform <- function(ds, ds_alt) {
  ds %>% 
    left_join(ds_alt, by=c("sex"="annot_sex")) %>% 
    mutate(text_sex=ifelse(is.na(alt_annot), "unknown", alt_annot)) %>%
    select(-sex, -alt_annot) %>%
    rename(expr_sex=sex_lab) %>%
    select(acc, text_sex, expr_sex, pred)
}
human2 <- addAltReform(human, human_alt)
mouse2 <- addAltReform(mouse, mouse_alt)
rat2 <- addAltReform(rat, rat_alt)

save(human2, mouse2, rat2, file="data/sample_level_sex.RData" )

# ----- agreement with metadata ----- #
h_filt <- human2 %>% 
  filter(!text_sex %in% c("unknown", "mixed")) 
sum(h_filt$text_sex==h_filt$expr_sex)/nrow(h_filt)

m_filt <- mouse2 %>% 
  filter(!text_sex %in% c("unknown", "mixed")) 
sum(m_filt$text_sex==m_filt$expr_sex)/nrow(m_filt)

r_filt <- rat2 %>% 
  filter(!text_sex %in% c("unknown", "mixed")) 
sum(r_filt$text_sex==r_filt$expr_sex)/nrow(r_filt)



# ------ set up study level ----- #
h3 <- h_map %>% full_join(human2, by=c("sample_acc"="acc")) %>% filter(!is.na(text_sex))
m3 <- m_map %>% full_join(mouse2, by=c("sample_acc"="acc"))  %>% filter(!is.na(text_sex))
r3 <- r_map %>% full_join(rat2, by=c("sample_acc"="acc")) %>% filter(!is.na(text_sex))


# ----- fraction missing ----- #

# sample level
sum(human2$text_sex=="unknown")/nrow(human2)
sum(mouse2$text_sex=="unknown")/nrow(mouse2)
sum(rat2$text_sex=="unknown")/nrow(rat2)

# ----- study level annot ----- #

# get study_type labels / counts
summarizeStudy <- function(ds){
 ds %>% 
  rename(gse=study_acc, gsm=sample_acc) %>%
  mutate(text_sex=ifelse(text_sex %in% c("unknown", "mixed"), NA, text_sex)) %>%
  rename(metadata=text_sex, exprsex=expr_sex) %>%
  gather(key="labeling_method", value="sex", -gse, -gsm) %>%
  group_by(gse, labeling_method) %>%
  dplyr::summarize(num_samples=n(),
            num_f=sum(sex=="female"),
            num_m=sum(sex=="male")) %>%
  mutate(study_type= case_when(
    (is.na(num_f) & is.na(num_m)) ~ "unlabeled",
    (!is.na(num_f) & !is.na(num_m) & num_f/num_samples > 0.8 & num_m > 0 ) ~ "mostly-female",
    (!is.na(num_f) & !is.na(num_m) & num_m/num_samples > 0.8 & num_f > 0 ) ~ "mostly-male",
    (!is.na(num_f) & !is.na(num_m) & num_f > 0 & num_m > 0 ) ~ "mixed",
    (!is.na(num_f) & num_f > 0 ) ~ "female-only",
    (!is.na(num_m) & num_m > 0 ) ~ "male-only")) %>%
  mutate(freq=1) %>% 
  ungroup(gse) %>%  
  mutate(study_type=factor(study_type, 
                           levels=c("female-only", "mostly-female", "mixed", "mostly-male", "male-only", "unlabeled")), 
         labeling_method=factor(labeling_method, 
                                levels=c("metadata", "exprsex")), 
         gse=as.factor(gse)) %>% 
  rename(study=gse) %>%
  rename(sex=study_type)
}

h_summ <- summarizeStudy(h3 %>% select(-pred))
m_summ <- summarizeStudy(m3 %>% select(-pred))
r_summ <- summarizeStudy(r3 %>% select(-pred))

h_summ %>% filter(labeling_method=="metadata") %>% group_by(sex) %>% count() 
length(unique(h_summ$study))

m_summ %>% filter(labeling_method=="metadata") %>% group_by(sex) %>% count() 
length(unique(m_summ$study))

r_summ %>% filter(labeling_method=="metadata") %>% group_by(sex) %>% count() 
length(unique(r_summ$study))

save(h_summ, m_summ, r_summ, file="data/summary_files.RData")

# ----- single sex ----- #

h_f <- h_summ %>% filter(labeling_method=="metadata" & sex=="female-only") %>% select(study)
h_f2 <- h3 %>% semi_join(h_f, by=c("study_acc"="study")) 
sum(h_f2$expr_sex==h_f2$text_sex)/nrow(h_f2)
h_m <- h_summ %>% filter(labeling_method=="metadata" & sex=="male-only") %>% select(study)
h_m2 <- h3 %>% semi_join(h_m, by=c("study_acc"="study")) 
sum(h_m2$expr_sex==h_m2$text_sex)/nrow(h_m2)


m_f <- m_summ %>% filter(labeling_method=="metadata" & sex=="female-only") %>% select(study)
m_f2 <- m3 %>% semi_join(m_f, by=c("study_acc"="study")) 
sum(m_f2$expr_sex==m_f2$text_sex)/nrow(m_f2)
m_m <- m_summ %>% filter(labeling_method=="metadata" & sex=="male-only") %>% select(study)
m_m2 <- m3 %>% semi_join(m_m, by=c("study_acc"="study")) 
sum(m_m2$expr_sex==m_m2$text_sex)/nrow(m_m2)


r_f <- r_summ %>% filter(labeling_method=="metadata" & sex=="female-only") %>% select(study)
r_f2 <- r3 %>% semi_join(r_f, by=c("study_acc"="study")) 
sum(r_f2$expr_sex==r_f2$text_sex)/nrow(r_f2)
r_m <- r_summ %>% filter(labeling_method=="metadata" & sex=="male-only") %>% select(study)
r_m2 <- r3 %>% semi_join(r_m, by=c("study_acc"="study")) 
sum(r_m2$expr_sex==r_m2$text_sex)/nrow(r_m2)
