# 01d_look_at_counts.R
# E Flynn
# 09/29/2020
# 
# Code for looking at the metadata sex breakdown. 
# This generates by-study and by-sample tables used for follow up analyses.
# We also plot the alluvial diagrams.
#
# TODO: 
# - when to exclude problem platforms?
# CONCERN: MOUSE RNA-seq shows a high fraction of female samples labeled ambiguous
#  - this indicates model misspecification!


require('tidyverse')
require('scales')  


#  --- set up color scheme --- #
library(RColorBrewer)
my.l <- brewer.pal(n = 8, name = 'Set2')
blues <- brewer.pal(9,name="Blues")
oranges <- brewer.pal(9,name="Oranges")
my.cols3 <- c (my.l[3],my.l[2], my.l[8])
my.cols4 <- c (my.l[3],my.l[4],my.l[2], my.l[8])
my.cols6 <- c (my.l[3], blues[4],my.l[4],oranges[4], my.l[2], my.l[8])

comb_metadata3 <- read_csv("data/data_old/01_sample_lists/combined_human_mouse_meta.csv", col_types= "cccccccdld")


# expression data QC + breakdown
# // TODO - expression data QC!
# - remove problem platforms


cutoff <- comb_metadata3 %>%
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


by_sample <- comb_metadata3 %>% 
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
ggsave("figures/paper_figs/figs1_sample_alluvial.png")  

# study level
study_flow_freq_counts <- by_study2 %>%
  mutate(organism=toupper(organism)) %>%
  mutate(data_type=ifelse(data_type=="rnaseq", "RNA-seq", "Microarray")) %>% 
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
  #geom_text(stat = "stratum", size = 4) +
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
  mutate(frac_unlabeled=unlabeled/num_samples,
         frac_female=female/num_samples,
         frac_male=male/num_samples) 


counts_sample_table %>%
  mutate(across(c(frac_unlabeled:frac_male), signif, 3)) %>%
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
  mutate(frac_unlabeled=unknown/num_studies,
         frac_male_only=`male only`/num_studies,
         frac_female_only=`female only`/num_studies,
         frac_mixed_sex=`mixed sex`/num_studies,
         frac_single_sex=(`male only`+`female only`)/num_studies)

# // TODO: format sig figs

counts_study_table %>% 
  mutate(across(contains("frac"), signif, 3)) %>%
  write_csv("tables/s1b_study_sex_counts.csv")


# are the differences significant? can't do this! 
# b/c RNA-seq metadata is lossy
counts_sample_table %>% 
  filter(label_type=="metadata") %>% 
  select(organism, data_type, unlabeled, num_samples) %>%
  mutate(num_present=num_samples-unlabeled) %>%
  select(-num_samples) 


# ---- Look at pooled/mixed sex samples ----- #

ggplot(comb_metadata3 %>% 
         filter(metadata_sex != "unknown") %>%
         mutate(metadata_sex=fct_recode(metadata_sex, "pooled"="mixed")) %>%
         rename(`metadata sex`=metadata_sex), 
       aes(x=p_male, col=`metadata sex`))+
  geom_density()+
  facet_grid(data_type~organism)+
  theme_bw()+
  scale_color_manual(values=my.cols3)+
  xlab("P(male)")
ggsave("figures/paper_figs/figs2_mixed_sex_dist.png")

frac_unlab <- comb_metadata3 %>%
  filter(!is.na(p_male) & metadata_sex!="unknown") %>%
  mutate(lab=case_when(
    p_male > 0.3 & p_male < 0.7 ~ "ambiguous",
    p_male <= 0.3 ~ "female",
    p_male >= 0.7 ~ "male")) %>%
  group_by(organism, data_type, metadata_sex,  lab) %>%
  count() %>%
  ungroup() %>%
  pivot_wider(names_from=lab, values_from=n) %>%
  group_by(organism, data_type, metadata_sex) %>%
  mutate(tot=ambiguous+female+male) %>%
  mutate(across(ambiguous:male, ~./tot)) %>%
  ungroup() %>%
  group_by(organism, data_type) %>%
  mutate(n=sum(tot), frac_tot=tot/n) 

frac_unlab %>%   
  mutate(across(c("ambiguous", "male", "female", "frac_tot"), signif, 3)) %>%
  rename("unclear"="ambiguous") %>%
  write_csv("tables/s1c_assignment_of_unclear.csv")
# CONCERN: MOUSE RNA-seq shows a high fraction of female sanokes labeled ambiguous
#  - this indicates model misspecification!


# ---- STOP ---- #
# --- Klinefelter's studies ---- #
klin_study <- comb_metadata %>% filter(str_detect(study_acc, "SRP117733"))
klin_study2 <- comb_metadata %>% filter(str_detect(study_acc, "GSE42331"))


klin_study %>%
  ggplot(aes(x=p_male))+geom_histogram()
klin_study %>%
  select(p_male, sample_acc) %>%
  arrange(p_male)

metadat <- read_csv("data/01_metadata/human_rnaseq_sample_metadata.csv", 
                    col_types=paste(rep("c", 11), collapse=""))
metadat2 <- read_csv("data/01_metadata/human_metadata.csv", 
                     col_types=paste(rep("c", 13), collapse=""))
klin_w_meta <- metadat %>% select(acc, title) %>% 
  inner_join(klin_study %>% select(sample_acc, p_male), 
             by=c("acc"="sample_acc")) %>%
  mutate(category=case_when(
    str_detect(title,"aKSL") ~ "Klin-like",
    str_detect(title, "aKS") ~ "Klin",
    str_detect(title, "pNorm") ~ "pre-pub norm",
    str_detect(title, "pKS") ~ "pre-pub Klin",
    str_detect(title, "aSCO") ~ "sertoli",
    str_detect(title, "aNorm") ~ "spermat"
  ))

klin_w_meta2 <- metadat2 %>% select(acc, title) %>% 
  inner_join(klin_study2 %>% select(sample_acc, p_male), 
             by=c("acc"="sample_acc")) %>%
  mutate(category=case_when(
    str_detect(title, "KS") ~ "KS",
    str_detect(title, "female") ~ "female",
    str_detect(title, "male") ~ "male"
  ))

ggplot(klin_w_meta, aes(y=category, x=p_male))+
  geom_point()
# aKS --> adult Klinefelter
# aKSL --> klinefelter like
# aSCO --> sertoli
# aNorm --> sperm
# pKS --> prepuberty Klin
# pNorm --> prepuberty norm

# concern: are we doing a bad job with prepuberty?

ggplot(klin_w_meta2, aes( x=p_male))+
  geom_histogram() +
  theme_bw() + 
  facet_grid(category~.)# this looks more reasonable but still not happy w
# ---------#

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




