# 05_drug_breakdown.R
# E Flynn
# Last updated 07/29/2020
#
# Look at the breakdown of drug data by sex. 
# TODO:
# - update so has the full set of studies! (including RNA-seq)
# - filter out cell line data
# - move ATC classes data + all other dependent data to this repo!
# - add in sample level data?
# - create supplemental tables with lists of drug studies

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


# --- construct the data frame --- #
# how we want the data formatted:
# organism | data_type | study | study_type | drug | drug_class
comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")
by_study <- read_csv("data/study_sex_lab.csv")


# this is missing RNA-seq!
drugbank_dat <- read_tsv("../labeling/geo2drug/data/02_labeled_data/drugbank_mapped_gse.txt")
# NOTE on drug breakdown data problems:
#   - 02_drug_gse_labeling.R (from geo2drug repo) - doesnt have RNA-seq 
#   - 02_map_to_drugbank.R (current repo) - sample level

study_db <- by_study %>% 
  left_join(drugbank_dat %>% select(gse, dbID,name, ATC), by=c("study_acc"="gse")) 

# --- plot sex breakdown of drug studies --- #
ggplot(study_db %>% 
         filter(!is.na(dbID)) %>%
         mutate(study_sex=factor(study_sex, levels=c("female only", "mostly female", 
                                                     "mixed sex", "mostly male", 
                                                     "male only", "unknown")),
                label_type=factor(label_type, levels=c("metadata", "expression"))), 
       aes(x=label_type))+
  geom_bar(aes(fill=study_sex))+
  facet_wrap(.~organism, scales="free")+
  ylab("Number of studies")+
  xlab("Label source")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=my.cols6)

ggsave("figures/paper_figs/supp_drug_sex_breakdown.png")

# table with counts for study breakdown
study_breakdown_drug <- study_db %>% 
  select(-name) %>%
  filter(label_type=="expression") %>%
  mutate(is_drug=(!is.na(dbID))) %>%
  group_by(organism, data_type, study_sex, is_drug) %>%
  count() %>%
  ungroup() %>%
  pivot_wider(id_cols=c(organism, data_type, study_sex), 
              names_from=is_drug, values_from=n, values_fill=list(n=0)) %>%
  rename(n_non_drug=`FALSE`, n_drug=`TRUE`)



# --- get sex breakdown at a sample level --- #

# doing single and multi-mapping separately for speed...
mult_map <- comb_metadata %>% 
  filter(str_detect(study_acc, ";")) %>%
  separate_rows(study_acc, sep=";") %>%
  left_join(drugbank_dat %>% select(gse, dbID), 
            by=c("study_acc"= "gse")) %>%
  select(-platform, -study_acc) %>%
  group_by(sample_acc) %>%
  mutate(dbID=paste(unique(dbID[!is.na(dbID)]), collapse=";")) %>%
  unique()

sing_map <- comb_metadata %>% 
  filter(!str_detect(study_acc, ";")) %>%
  left_join(drugbank_dat %>% select(gse, dbID), 
            by=c("study_acc"= "gse")) %>%
  select(-platform, -study_acc) 

drug_sample_dat <- mult_map %>% 
  ungroup() %>%
  bind_rows(sing_map) %>%
  mutate(is_drug=(!is.na(dbID) & dbID != ""))
  
# plot drug/non-drug vs sample sex
ggplot(data.frame(drug_sample_dat) %>%
         mutate(expr_sex=ifelse(is.na(expr_sex), "unknown", expr_sex)),
       aes(x=is_drug))+
  geom_bar(aes(fill=expr_sex))+
  facet_grid(data_type~organism, scales="free")+
  ylab("Number of samples")+
  xlab("is drug?")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=my.cols3)+
  scale_y_continuous(labels = comma) 


ggsave("figures/paper_figs/supp_drug_sample_breakdown_v0.png")

drug_sample_dat %>% 
  filter(data_type=="microarray") %>%
  group_by(organism, is_drug, expr_sex) %>%
  count()


# ---- Make figure with ATC breakdown
# add ATC descriptions
atc_names <- read.delim("../labeling/geo2drug/data/deprecated/ref_data/atc_classes.txt",
                        header=FALSE, stringsAsFactors = FALSE)
colnames(atc_names) <- c("class", "descript")

study_db2 <- study_db %>% 
  filter(!is.na(dbID)) %>%
  mutate(class=substr(ATC, 1, 1)) %>%
  left_join(atc_names, by=c("class"="class")) %>%
  filter(label_type!="metadata") %>%
  mutate(study_sex=case_when(
    study_sex=="mostly female" ~ "mixed sex",
    study_sex=="mostly male" ~ "mixed sex", 
    TRUE ~ study_sex)) %>%
  filter(!is.na(descript)) %>%
  mutate(ATC_class=paste(class, descript, sep=" - ")) %>%
  select(-label_type, -descript)

study_db2 %>% write_csv("data/drug_studies.csv")

# plot ATC breakdown
ggplot(study_db2 %>% filter(organism=="human"), 
       aes(x=class, fill=ATC_class))+
  geom_histogram(stat="count")+
  facet_grid(study_sex ~ .)+  
  ylab("Number of studies")+
  xlab("ATC class")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())
ggsave("figures/paper_figs/fig5_drug_breakdown_human.pdf")

ggplot(study_db2 %>% filter(organism=="mouse"), 
       aes(x=class, fill=ATC_class))+
  geom_histogram(stat="count")+
  facet_grid(study_sex ~ .)+  
  ylab("Number of studies")+
  xlab("ATC class")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())
ggsave("figures/paper_figs/fig5_drug_breakdown_mouse.pdf")


# --- Run statistical tests for enrichment --- #
count_per_class_df <- study_db2 %>% 
  select(organism, class, study_sex) %>%
  group_by(organism, study_sex) %>%
  mutate(tot=n()) %>%
  group_by(organism, class, study_sex) %>% 
  mutate(num_in_grp=n(),
         num_not_in_grp=tot-num_in_grp) %>%
  unique() %>%
  ungroup()  %>%
  arrange(organism, class, study_sex) 

# COMPARISONS:
# - single vs mixed sex
# - male vs female only

# create lists of data frames, grouped by organism and class  
class_mixed_vs_single <- count_per_class_df %>%
  mutate(study_sex=case_when(
    study_sex=="male only" ~ "single sex",
    study_sex=="female only" ~ "single sex",
    study_sex=="mixed sex" ~ "mixed sex"
  )) %>%
  group_by(organism, class, study_sex) %>%
  summarize(num_in_grp=sum(num_in_grp),
            num_not_in_grp=sum(num_not_in_grp)) %>%
  ungroup() %>%
  group_split(organism, class) 
  
class_f_vs_m <- count_per_class_df %>%
  filter(study_sex!="mixed sex") %>%
  group_split(organism, class) 
  
class_names <- lapply(class_f_vs_m, function(x) 
  paste(x$organism[[1]], x$class[[1]], sep="_"))

# function for running a chisq test on a data frame
runCHisq <- function(grp_data, class_names) {
  res_dat <- lapply(grp_data, function(x) { 
  my_df <- data.frame(x[,c("num_in_grp", "num_not_in_grp")])
  too_small <- (min(my_df) <= 5)
  if (too_small) {
    return(NA)
  }
  rownames(my_df) <- x$study_sex
  res <- chisq.test(my_df)
  pval <- res$p.value
  resid <- res$residuals[,"num_in_grp"]
  return(list("pval"=pval, "resid1"=resid[1], "resid2"=resid[2] ))
  })
  names(res_dat) <- class_names
  res_df <- data.frame(do.call(rbind,res_dat[!is.na(res_mf)]))
  res_df$grp <- rownames(res_df)
  return(res_df)
}

# run chisq tests
fm_chisq <- runCHisq(class_f_vs_m, class_names)
ms_chisq <- runCHisq(class_mixed_vs_single, class_names)

# filter for pval cutoff 
(p_cut <- 0.05/(nrow(ms_chisq)+nrow(fm_chisq)) ) # p < 0.001
ms_chisq %>% filter(pval < p_cut) %>%
  rename("resid_mixed_sex"=resid1, "resid_single_sex"=resid2)

fm_chisq %>% filter(pval < p_cut) %>%
  rename("resid_f"=resid1, "resid_m"=resid2)

# ------ COMPARE TO WANG ET AL ------ #
drug_pert_manual <- read_csv("../drug_expression/drug_labeling/data/single_drug_perturbations-v1.0.csv")
drug_pert_auto <- read_csv("../drug_expression/drug_labeling/data/single_drug_perturbations-p1.0.csv")
length(unique(drug_pert_manual$geo_id)) # 337
length(unique(drug_pert_auto$geo_id)) # 378



# add in sex labels
manual_sl <- drug_pert_manual %>% 
  filter(organism != "rat") %>%
  select(drug_name, drugbank_id, geo_id, organism) %>%
  unique() %>%
  rename(study_acc=geo_id) %>%
  left_join(by_study  %>% 
              filter(label_type=="expression")%>% 
              select(-organism, -data_type, -label_type), by=c("study_acc")) 
for_plot <- manual_sl %>% 
  select(study_acc, study_sex, organism) %>% 
  unique() %>%
  mutate(study_sex=ifelse(is.na(study_sex), "unknown", study_sex)) %>%
  mutate(study_sex=factor(study_sex, levels=c("female only", "mostly female", 
                                              "mixed sex", "mostly male", 
                                              "male only", "unknown")))

ggplot(for_plot, aes(x=organism, fill=study_sex)) + 
  geom_bar()+
  ylab("Number of studies")+
  xlab("")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=my.cols6)
ggsave("figures/paper_figs/wang_manual_breakdown.png")

# add in ATC
drug_info <- read.delim("../labeling/geo2drug/data/00_db_data/drugbank_parsed.txt", 
                             stringsAsFactors = FALSE) %>%
  select(name, dbID, ATC)

manual_sl_class <- manual_sl %>% 
  left_join(drug_info, by=c("drugbank_id"="dbID")) %>%
  filter(!is.na(drugbank_id)) %>%
  mutate(class=substr(ATC, 1, 1)) %>%
  left_join(atc_names, by=c("class"="class")) %>%
  mutate(study_sex=case_when(
    is.na(study_sex) ~ "unknown",
    study_sex=="mostly female" ~ "mixed sex",
    study_sex=="mostly male" ~ "mixed sex", 
    TRUE ~ study_sex)) %>%
  filter(!is.na(descript)) %>%
  mutate(ATC_class=paste(class, descript, sep=" - ")) %>%
  select( -descript)

ggplot(manual_sl_class, 
       aes(x=class, fill=ATC_class))+
  geom_histogram(stat="count")+
  facet_grid(study_sex ~ organism)+  
  ylab("Number of studies")+
  xlab("ATC class")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())
ggsave("figures/paper_figs/wang_manual_atc_class.png")
# // TODO: look at automated data?


# ----- by drug breakdown ---- #
study_types <- study_db2 %>%  
  group_by(organism, data_type, name, study_sex) %>% 
  count() %>%
  ungroup() %>%
  pivot_wider(id_cols=c(organism, data_type, name), 
              names_from="study_sex", values_from="n",
              values_fill=list(n=0)) %>%
  arrange(organism, data_type, name) %>%
  mutate(study_types=case_when(
    `mixed sex` !=0 & `male only` == 0 & `female only`==0 ~ "mixed sex",
    `mixed sex`==0 & `male only` != 0 & `female only`!=0 ~ "single sex - both",
    `mixed sex`!=0 & (`male only` != 0 | `female only`!=0 )~ "both",
    `mixed sex`==0 & `male only` == 0 & `female only`!=0 ~ "single sex - f",
    `mixed sex`==0 & `male only` != 0 & `female only`==0 ~ "single sex - m",
    TRUE ~ "other"
    ))

# number of drugs per study type
count_drugs <- study_types %>% 
  group_by(organism, data_type, study_types) %>% 
  count() %>% 
  ungroup() %>%
  group_by(organism, data_type) %>% 
  mutate(sum=sum(n)) %>%
  ungroup() %>%
  mutate(frac=n/sum)

# ------ LOOK AT CANCER DRUGS ------ #
cancer_drugs <- study_db2 %>% filter(class=="L", organism=="human") 

sex_breakdown_cancer_drugs <- cancer_drugs %>% 
  group_by(organism, data_type, dbID, name, study_sex) %>%
  mutate(study_acc=paste(study_acc, collapse = ";")) %>%
  select(-num_tot, -num_f, -num_m, -num_present) %>%
  unique() %>%
  ungroup() %>%
  pivot_wider(id_cols=c(organism, data_type, name, dbID), 
              names_from="study_sex", values_from="study_acc") %>%
  arrange(organism, data_type, name) 
sex_breakdown_cancer_drugs2 <- sex_breakdown_cancer_drugs %>%
  group_by(organism, data_type, name) %>%
  mutate(across(c(`mixed sex`, `male only`, `female only`), ~ifelse(is.na(.), 0, length(str_split(., ";")[[1]])))) %>%
  mutate(study_types=case_when(
    `mixed sex` !=0 & `male only` == 0 & `female only`==0 ~ "mixed sex",
    `mixed sex`==0 & `male only` != 0 & `female only`!=0 ~ "single sex - both",
    `mixed sex`!=0 & (`male only` != 0 | `female only`!=0 )~ "both",
    `mixed sex`==0 & `male only` == 0 & `female only`!=0 ~ "single sex - f",
    `mixed sex`==0 & `male only` != 0 & `female only`==0 ~ "single sex - m",
    TRUE ~ "other"
  )) %>%
  ungroup()

sex_breakdown_cancer_drugs2 %>% 
  group_by(study_types) %>% 
  count()

sex_breakdown_cancer_drugs2 %>% 
  mutate(tot=`male only`+`female only`+`mixed sex`) %>% 
  arrange(desc(tot))
cancer_drugs %>% 
  filter(organism=="human", data_type =="microarray") %>% 
  select(study_acc, study_sex)  %>% unique() %>%
  group_by(study_sex) %>% count()

# gotten count data, want more of a narrative
# REALLY NEED TO FILTER OUT CELL LINES

# ------ LOOK AT NEURO DRUGS ------- #
neuro_drugs <- study_db2 %>% filter(class=="N")

sex_breakdown_neuro_drugs <- neuro_drugs %>% 
  group_by(organism, data_type, dbID, name, study_sex) %>%
  mutate(study_acc=paste(study_acc, collapse = ";")) %>%
  select(-num_tot, -num_f, -num_m, -num_present) %>%
  unique() %>%
  ungroup() %>%
  pivot_wider(id_cols=c(organism, data_type, name, dbID), 
              names_from="study_sex", values_from="study_acc") %>%
  arrange(organism, data_type, name) 

#  mouse nervous system drugs and associated studies
sex_breakdown_neuro_drugs %>% write_csv("tables/sex_breakdown_neuro_drugs_by_drug.csv")

sex_breakdown_neuro_drugs2 <- sex_breakdown_neuro_drugs %>%
  mutate(across(c(`mixed sex`, `male only`, `female only`), ~ifelse(is.na(.), 0, 1))) %>%
  mutate(study_types=case_when(
    `mixed sex` !=0 & `male only` == 0 & `female only`==0 ~ "mixed sex",
    `mixed sex`==0 & `male only` != 0 & `female only`!=0 ~ "single sex - both",
    `mixed sex`!=0 & (`male only` != 0 | `female only`!=0 )~ "both",
    `mixed sex`==0 & `male only` == 0 & `female only`!=0 ~ "single sex - f",
    `mixed sex`==0 & `male only` != 0 & `female only`==0 ~ "single sex - m",
    TRUE ~ "other"
  ))

neuro_drugs_count <- sex_breakdown_neuro_drugs2 %>% 
  group_by(organism, data_type, study_types) %>% 
  count() %>% 
  ungroup() %>%
  group_by(organism, data_type) %>% 
  mutate(sum=sum(n)) %>%
  ungroup() %>%
  mutate(frac=n/sum)

neuro_drugs_count
m_only_mouse <-sex_breakdown_neuro_drugs2 %>% 
  filter(organism=="mouse" & study_types == "single sex - m") 
