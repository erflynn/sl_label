# 05a_drug_breakdown.R
# E Flynn
# Last updated 07/29/2020
#
# Look at the breakdown of drug data by sex. 
# TODO:
# - filter out cell line data
# - move ATC classes data + all other dependent data to this repo!
# - add in sample level data?
# - create supplemental tables with lists of drug studies

require('tidyverse')
require('scales')
require('ggrepel')


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
comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", col_types="cccccccdld")
by_study <- read_csv("data/study_sex_lab.csv")

# load drugbank data 
# each row is a study/drug not a study!
drugbank_study_dat <- read_csv("data/02_labeled_data/study_to_drug.csv")
drugbank_study_level <- drugbank_study_dat %>%
  group_by(organism, data_type, study_acc) %>%
  summarise(across(dbID:ATC, ~paste(unique(.), collapse=";"))) %>%
  ungroup()


study_db <- by_study %>% 
  left_join(drugbank_study_level %>% 
              select(study_acc, dbID, name, ATC), by=c("study_acc")) 

# --- plot sex breakdown of drug studies --- #
ggplot(study_db %>% 
         filter(!is.na(dbID)) %>%
         mutate(study_sex=factor(study_sex, levels=c("female only", "mostly female", 
                                                     "mixed sex", "mostly male", 
                                                     "male only", "unknown")),
                label_type=factor(label_type, levels=c("metadata", "expression"))), 
       aes(x=label_type))+
  geom_bar(aes(fill=study_sex))+
  facet_grid(data_type~organism, scales="free")+
  ylab("Number of studies")+
  xlab("Label source")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=my.cols6)

ggsave("figures/paper_figs/supp_drug_sex_breakdown.png")

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
  left_join(drugbank_study_dat %>% select(study_acc, dbID), 
            by=c("study_acc")) %>%
  select(-platform, -study_acc) %>%
  group_by(sample_acc) %>%
  mutate(dbID=paste(unique(dbID[!is.na(dbID)]), collapse=";")) %>%
  unique()

sing_map <- comb_metadata %>% 
  filter(!str_detect(study_acc, ";")) %>%
  left_join(drugbank_study_dat %>% 
              select(study_acc, dbID), 
            by=c("study_acc")) %>%
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

study_db2 <- by_study %>% 
  left_join(drugbank_study_dat %>%  # join to study/drug
              select(study_acc, dbID, name, ATC), by=c("study_acc"))%>%
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
  select(-label_type, -descript) %>%
  select(organism, data_type, study_acc, everything())  %>%
  arrange(organism, data_type, study_acc)

study_db2 %>% write_csv("data/drug_studies.csv") # study-drug is a row

# plot ATC breakdown
ggplot(study_db2 %>% filter(organism=="human" & study_sex!="unknown"), 
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

ggplot(study_db2 %>% filter(organism=="mouse" & study_sex!="unknown"), 
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
#counts_per_class_df <- manual_sl_class %>%
#count_per_class_df <- hc_drugs_dat %>% 
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
  filter(!study_sex %in% c("unknown")) %>%
  mutate(study_sex=case_when(
    study_sex=="male only" ~ "single sex",
    study_sex=="female only" ~ "single sex",
    study_sex=="mixed sex" ~ "mixed sex"
  )) %>%
  group_by(organism, class, study_sex) %>%
  summarise(num_in_grp=sum(num_in_grp),
            num_not_in_grp=sum(num_not_in_grp)) %>%
  ungroup() %>%
  group_split(organism, class) 
  
class_f_vs_m <- count_per_class_df %>%
  filter(!study_sex %in% c("unknown")) %>%
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
  print(is.na(res_dat))
  print(class_names[!is.na(res_dat)])
  
  #names(res_dat) <- class_names
  res_df <- data.frame(do.call(rbind, res_dat[!is.na(res_dat)]))
  
  res_df$grp <- class_names[!is.na(res_dat)]
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

# ------ COMPARE TO HC STUDIES ------ #
hc_drugs <- read_csv("data/hc_drug_labels.csv")
hc_drugs %>%
  left_join(by_study  %>% 
              filter(label_type=="expression")%>% 
              select(-organism, -data_type, -label_type), by=c("study_acc")) %>%
  select(study_acc, study_sex, organism) %>% 
  unique() %>%
  mutate(study_sex=ifelse(is.na(study_sex), "unknown", study_sex)) %>%
  mutate(study_sex=factor(study_sex, levels=c("female only", "mostly female", 
                                              "mixed sex", "mostly male", 
                                              "male only", "unknown"))) %>%
  ggplot(aes(x=organism, fill=study_sex)) + 
  geom_bar()+
  ylab("Number of studies")+
  xlab("")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=my.cols6)

hc_drugs_dat <- hc_drugs %>%
  left_join(by_study  %>% 
              filter(label_type=="expression")%>% 
              select(-organism, -data_type, -label_type), by=c("study_acc")) %>%
  mutate(study_sex=ifelse(is.na(study_sex), "unknown", study_sex)) %>%
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
hc_drugs_dat %>%  ggplot(aes(x=class, fill=ATC_class))+
  geom_histogram(stat="count")+
  facet_grid(study_sex ~ organism)+  
  ylab("Number of studies")+
  xlab("ATC class")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

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

length(unique(manual_sl$study_acc))
length(unique(manual_sl$drugbank_id))

# what is the overlap?
manual_geo_studies <- unique(manual_sl$study_acc) # 299
auto_geo_studies <- drug_pert_auto %>% filter(organism != "rat") %>%
  distinct(geo_id) %>% pull() # 308

our_drug_studies <- drugbank_study_level %>% 
  filter(data_type=="microarray" & 
           str_detect(study_acc, "GSE")) %>%
  distinct(study_acc) %>% pull() # 4461

length(intersect(manual_geo_studies, our_drug_studies)) # 195
length(intersect(auto_geo_studies, our_drug_studies)) # 122

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

# ------ plot sex diff in drugs ------ #
all_no_cl <- study_db2 %>% 
  inner_join(stud_tiss2 %>% 
               select(study_acc, tissue, cancer, cell_line, tot )) %>%
  filter(tissue+cancer > 0.5)

sex_breakdown_drugs <- all_no_cl %>% 
  select(-tissue, -cell_line, -tot, -cancer) %>%
  group_by(organism,  dbID, class, name, study_sex) %>%
  mutate(study_acc=paste(study_acc, collapse = ";")) %>%
  select(-num_tot, -num_f, -num_m, -num_present, -data_type) %>%
  unique() %>%
  ungroup() %>%
  pivot_wider(id_cols=c(organism, name, dbID, class), 
              names_from="study_sex", values_from="study_acc") %>%
  arrange(organism, name) 

sex_breakdown_drugs2 <- sex_breakdown_drugs %>%
  group_by(organism, name, class) %>%
  mutate(across(c(`mixed sex`, `male only`, `female only`), 
                ~ifelse(is.na(.), 0, length(str_split(., ";")[[1]])))) %>%
  mutate(num_studies=`mixed sex`+`female only`+`male only`) %>%
  select(-unknown) %>%
  mutate(across(c(`mixed sex`, `female only`, `male only`), ~./num_studies))
sex_breakdown_drugs2 %>% write_csv("data/by_drug_study_sex_fraction.csv")
sex_breakdown_drugs2 %>%
  mutate(drug_name=case_when(
    `female only` > 0.66 & `male only` < 0.33 & num_studies >= 3 ~name,
    `female only` < 0.33 & `male only` > 0.66 & num_studies >= 3 ~ name,
    TRUE ~ ""))  %>%
  rename("number of studies per drug"=num_studies) %>%
  ggplot(aes(x=`female only`, y=`male only`, col=class))+
  geom_point(alpha=0.5, aes(size=`number of studies per drug`), 
             position=position_jitter())+
  geom_label_repel(size=3, fill="white", aes(label=drug_name))+
  theme_bw()+
  xlab("proportion of female only studies")+
  ylab("proportion of male only studies")+
  facet_grid(.~organism)
ggsave("figures/paper_figs/sex_breakdown_drugs.png")

# ------ LOOK AT CANCER DRUGS ------ #
cancer_drugs <- study_db2 %>% filter(class=="L") 


sample_source <- read_csv("data/sample_source_type.csv") %>%
  select(acc, source_type)

stud_tiss <- comb_metadata %>% select(sample_acc, study_acc) %>%
  inner_join(sample_source, by=c("sample_acc"="acc")) %>%
  separate_rows(study_acc, sep=";") %>%
  mutate(source_type=fct_collapse(source_type,
    "cell_line"=c("unnamed_cl", "named_cl"),
    "other"=c("xenograft", "other")
  )) %>%
  group_by(study_acc, source_type) %>% 
  count() 

stud_tiss2 <- stud_tiss %>%
  pivot_wider(names_from=source_type, values_from=n, values_fill=0) %>%
  mutate(tot=tissue+primary_cells+cell_line+cancer+other+stem_cell) %>%
  mutate(across(tissue:stem_cell, ~./tot))

# ---- GRAB MIXED SEX PRIMARY AND STEM CELLS ---- #
by_study %>% filter(label_type=="expression") %>%
                          select(study_acc, organism, data_type, study_sex, num_f, num_m) %>% 
  filter(study_sex=="mixed sex" ) %>%
  semi_join(stud_tiss2 %>% filter(primary_cells==1 & tot >= 8)) %>%
  arrange(desc(num_f+num_m))


by_study %>% filter(label_type=="expression") %>%
  select(study_acc, organism, data_type, study_sex, num_f, num_m) %>% 
  filter(study_sex=="mixed sex" ) %>%
  semi_join(stud_tiss2 %>% filter(stem_cell==1 & tot >= 8)) %>%
  arrange(desc(num_f+num_m))
# ----------------------------------------------------- #
  
stud_tiss %>% group_by(source_type) %>% summarise(sum(n))
cc_no_cl <- cancer_drugs %>% 
  inner_join(stud_tiss2 %>% 
               select(study_acc, tissue, cancer, cell_line, tot )) %>%
  filter(tissue+cancer > 0.5)

cc_no_cl %>% filter(name=="Fulvestrant", organism=="human", data_type=="microarray") %>% 
                    select(-dbID, -ATC, -class, -ATC_class)
# GSE27444 - MCF7
# GSE33658 - patients
# GSE48905 - patients

cc_no_cl %>% filter(name=="Docetaxel", organism=="human", data_type=="microarray") %>% 
  select(-dbID, -ATC, -class, -ATC_class)
# all of these are looking at breast cancer

cc_no_cl %>% filter(name=="Decitabine", organism=="human") %>% 
  select(-dbID, -ATC, -class, -ATC_class)
# GSE19610 (m) - MDS stem cells
# GSE44857 (m) - AML xenografts 
# GSE55410 - ovarian
# SRP114518 - ovarian
# SRP046233 (m) - 14 CMML patients
# remaining - cell lines :/ 

sex_breakdown_cancer_drugs <- cc_no_cl %>% 
  select(-tissue, -cell_line, -tot, -cancer) %>%
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

sex_breakdown_cancer_drugs2 %>% arrange(desc(`female only`))

cancer_drugs3 <- sex_breakdown_cancer_drugs2 %>%
  select(-unknown, -study_types) %>%
  #filter(data_type=="microarray", organism=="mouse") %>%
  mutate(num_studies=`mixed sex`+`female only`+`male only`) %>%
  mutate(across(`mixed sex`:`male only`, ~./num_studies)) %>%
  mutate(drug_name=case_when(
    `female only` > 0.5 & `male only` < 0.5 & num_studies >= 2 ~name,
    `female only` < 0.5 & `male only` > 0.5 & num_studies >= 2 ~ name,
    TRUE ~ "")) 
cancer_drugs3 %>%
  ggplot(aes(x=`female only`, y=`male only`, 
              label=drug_name))+
    geom_point(alpha=0.5, aes(size=num_studies))+
  geom_label_repel(size=2, fill="white")+
  theme_bw()+
  facet_grid(organism ~ data_type)

cancer_type <- read_tsv("tables/cancer_type_annot.txt")
cancer_types2 <- cancer_type %>%
  mutate(cancer_type2=case_when(
    cancer_type=="breast" ~ "breast cancer",
    cancer_type=="prostate" ~ "prostate cancer",
    cancer_type=="ovarian" ~ "ovarian cancer",

    str_detect(cancer_type, "breast") &
      str_detect(cancer_type, "ovarian") & 
      str_detect(cancer_type, "prostate") ~ "includes breast,ovarian,prostate",
    str_detect(cancer_type, "ovarian") & 
      str_detect(cancer_type, "prostate") ~ "includes ovarian,prostate",
    str_detect(cancer_type, "breast") & 
      str_detect(cancer_type, "prostate") ~ "includes breast,prostate",
    str_detect(cancer_type, "breast") & 
      str_detect(cancer_type, "ovarian") ~ "includes breast,ovarian",
    str_detect(cancer_type, "breast") ~ "includes breast",
    
    str_detect(cancer_type, "prostate") ~ "includes prostate",
    str_detect(cancer_type, "ovarian") ~ "includes ovarian",
    TRUE ~ "other"
    ))

cancer_types2 %>% 
  group_by(cancer_type2) %>% count() %>% 
  arrange(desc(n))

cancer_drugs4 <- cancer_drugs3 %>% 
  left_join(cancer_types2 %>% select(-cancer_type)) %>%
  rename(cancer_type=cancer_type2) %>%
  mutate(cancer_type=ifelse(is.na(cancer_type), "not labeled", cancer_type))
cancer_drugs4 %>%  filter(drug_name!="") %>%
  arrange(organism, data_type, desc(`female only`), desc(`num_studies`)) %>%
  select(-`mixed sex`, -drug_name) %>%
  write_csv("tables/sex_breakdown_cancer_drugs.csv")

cancer_drugs4 %>%
  rename(`number of studies`=`num_studies`,
         `cancer type`=`cancer_type`) %>%
  ggplot(aes(x=`female only`, y=`male only`, 
             label=drug_name))+
  geom_point(alpha=0.7, aes(size=`number of studies`, col=`cancer type`))+
  geom_label_repel(size=2, fill="white", aes(col=`cancer type`))+
  theme_bw()+
  scale_color_manual(values=c("purple", "pink", "red", "green", 
                              "orange", "gray","black", "blue"))+
  facet_grid(organism ~ data_type)

cancer_drugs4 %>%   
  filter(drug_name!="") %>%
  arrange(organism, data_type, desc(`female only`), desc(`num_studies`)) %>%
  select(-`mixed sex`, -drug_name) %>%
  View()

# ------ LOOK AT NEURO DRUGS ------- #
neuro_drugs <- study_db2 %>% filter(class=="N")

neuro_drugs %>%
  group_by(organism, data_type, dbID, name, study_sex) %>%
  count() %>%
  pivot_wider(names_from="study_sex", values_from="n", values_fill=0) %>%
  select(-unknown) %>%
  filter(organism=="mouse") %>%
  mutate(num_studies=`female only`+`male only`+`mixed sex`) %>%
  mutate(across(`female only`:`mixed sex`, ~./num_studies)) %>%
  arrange(data_type, desc(num_studies), desc(`male only`), desc(`female only`)) %>%
  write_csv("tables/mouse_neuro_drugs.csv")

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



