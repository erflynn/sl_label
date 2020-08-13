# 01d_check_study_acc.R
# E Flynn
# 8/13/2020
#
# Code for looking at study-level accuracy
# Looks at study-level accuracy for single and mixed sex studies
# in the extended test set
#
# Makes the following figures:
#  - study-level accuracy/fraction labeled boxplot
#  - confusion matrix of metadata vs predicted sex
#
# TODO 
# - expand confusion matrix to consider include "mostly"

require('tidyverse')

sample_breakdown <- function(ds){
  ds %>%
    mutate(unlab=(p_male < 0.7 & p_male > 0.3 & !is.na(p_male))) %>%
    group_by(study_acc) %>%
    summarize(
      num_samples=n(),
      num_f=sum(metadata_sex=="female" & !is.na(p_male)),
      num_m=sum(metadata_sex=="male" & !is.na(p_male)),
      num_missing=sum(is.na(p_male)),
      num_unlab=sum(unlab, na.rm=TRUE),
      num_correct=sum(metadata_sex==expr_sex & !unlab, na.rm=TRUE)) %>%
    mutate(num_present=num_samples-num_missing,
           accuracy=num_correct/(num_present-num_unlab)) %>%
    mutate(frac_labeled=(num_present-num_unlab)/num_present) 
  
}



study_level <- function(my_organism, my_data_type){
  extended_test2 <- read_csv(sprintf("data/%s_%s_extended_test.csv", 
                                     my_organism, my_data_type))
  
  test_studies <- extended_test2 %>% 
    distinct(study_acc) %>% 
    separate_rows(study_acc, sep=";") %>% 
    unique() %>% pull(study_acc)
  by_study_sm <- by_study %>% 
    filter(label_type=="metadata" & !study_sex=="unknown" &
             organism==my_organism & data_type==my_data_type &
             study_acc %in% test_studies) 
  
 
  
  ms <- by_study_sm %>% 
    filter(study_sex == "mixed sex" & num_tot >= 10) 
  ss_f <- by_study_sm %>% 
    filter(num_tot==num_f & num_tot >= 8) 
  ss_m <- by_study_sm %>% 
    filter(num_tot==num_m & num_tot >= 8) 
  
  
  ms2 <- addSamples(ms, extended_test2)
  ss_f2 <- addSamples(ss_f, extended_test2)
  ss_m2 <- addSamples(ss_m, extended_test2)
  ms3 <- ms2 %>% sample_breakdown()
  ssf3 <- ss_f2 %>% sample_breakdown() 
  ssm3 <- ss_m2 %>% sample_breakdown() 
  my_dat <- bind_rows(ms3 %>% mutate(ds="mixed_sex"),
                      ssf3 %>% mutate(ds="single_sex_f"),
                      ssm3 %>% mutate(ds="single_sex_m")) %>%
    mutate(organism=my_organism, data_type=my_data_type)
  
  # add in the study sex
  
  return(my_dat)
}
# what is the accuracy/frac_labeled broken down by study
hr_study <- study_level("human", "rnaseq")
mr_study <- study_level("mouse", "rnaseq")
hm_study <- study_level("human", "microarray")
mm_study <- study_level("mouse", "microarray")

study_dat <- 
  bind_rows(hr_study, mr_study, mm_study, hm_study) 
study_dat2 <- study_dat %>%
  select(organism, data_type, ds, num_present, accuracy, frac_labeled, study_sex) 

study_dat2 %>%
  pivot_longer(cols=c("accuracy", "frac_labeled"), 
               names_to="metric") %>%
  mutate(ds=fct_recode(ds, 
                       "mixed sex" = "mixed_sex",
                       "female only"="single_sex_f",
                       "male only"="single_sex_m"),
         metric=fct_recode(metric,
            "fraction labeled"="frac_labeled")) %>%
  mutate(num_present=ifelse(num_present > 200, 200, num_present)) %>%
  rename("number of samples in study"="num_present") %>%
  unite(grp, c(metric, ds), remove=FALSE) %>%
  ggplot(aes(x=ds, y=value, group=grp, col=metric))+
  geom_boxplot()+
  geom_point(aes(size=`number of samples in study`), alpha=0.3, 
             position=position_jitterdodge())+
  theme_bw()+
  ylab("")+
  xlab("")+
  facet_grid(data_type ~ organism)+
ggsave("figures/paper_figs/study_acc_at_cutoff.png")  
  

expr_study_sex <- by_study %>% 
  filter(label_type=="expression") %>% 
  distinct(study_acc, study_sex) 
study_dat3 <- study_dat %>% 
  select(-study_sex) %>%
  separate_rows(study_acc, sep=";") %>%
  left_join(expr_study_sex)

# what fraction of study sex is accurate
study_dat4 <- study_dat3 %>% 
  select(organism, data_type, ds, study_sex) %>%
  group_by(organism, data_type, ds) %>%
  mutate(tot=n()) %>%
  ungroup() %>%
  mutate(ds=fct_recode(ds,
                       "mixed sex"="mixed_sex",
                       "female only"="single_sex_f",
                       "male only"="single_sex_m"
  )) %>%
    group_by(organism, data_type, ds, study_sex, tot) %>% 
  count() %>%
  pivot_wider(id_cols=c("organism", "data_type", 
                        "ds", "tot"),
              values_from=n, 
              names_from=study_sex, 
              values_fill=0) %>%
  pivot_longer(c("female only", 
                 "mostly female", 
                 "mixed sex", 
                 "mostly male",
                 "male only", 
                 "unknown"),
               values_to="n", 
               names_to="study_sex") %>%
  mutate(ds=fct_relevel(ds,
                        c("female only", 
                          "mixed sex", 
                          "male only")),
         study_sex=fct_relevel(study_sex, 
                               c("female only", 
                                 "mostly female", 
                                 "mixed sex", 
                                 "mostly male",
                                 "male only", 
                                 "unknown"))) %>%
  mutate(`fraction of dataset`=n/tot) 

# make a confusion matrix   
ggplot(study_dat4, aes(x=ds, y=study_sex))+
  geom_tile(aes(fill = `fraction of dataset`), colour = "white") +
  geom_text(aes(label = n)) +
  scale_fill_gradient(low="light blue", high="blue") +
  theme_bw() +
  ylab("predicted study sex")+
  xlab("metadata study sex")+
  facet_grid(organism ~ data_type)
ggsave("figures/paper_figs/study_confusion_matrix.png")


study_dat4 %>% 
  ungroup() %>% select(-tot, -n) %>%
  mutate(across(where(is.factor), as.character)) %>%
  filter(study_sex==ds) %>%
  select(-ds) %>%
  pivot_wider(names_from=study_sex, values_from=`fraction of dataset`)
