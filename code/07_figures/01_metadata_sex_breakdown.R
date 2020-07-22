# 01_metadata_sex_breakdown.R
# E Flynn
# 7/20/2020
# code for looking at the metadata sex breakdown 
# + also producing the table with the missingness -AND- overlap

require('tidyverse')
require('data.table')
require('scales')  


# # ----- 1. OVERLAP ----- #
# 
# # S0a/b
# m_overlap_ids <- intersect(m_r$acc, m_m$acc) 
# m_microarray <- setdiff(m_m$acc, m_r$acc)
# m_rnaseq_ds_only <- setdiff(m_r$acc, m_m$acc) 
# h_overlap_ids <- intersect(h_r$acc, h_m$acc) # 99611 overlapping IDs
# h_microarray <- setdiff(h_m$acc, h_r$acc) # 330508
# h_rnaseq_ds_only <- setdiff(h_r$acc, h_m$acc) # 130178
# 
# 
# require('VennDiagram')
# grid.newpage()
# draw.pairwise.venn(
#   area1=length(h_r$acc),
#   area2=length(h_m$acc),
#   cross.area=length(h_overlap_ids),
#   category=c(sprintf("RNA-seq\n(%s)", length(h_r$acc)), 
#                    sprintf("compendia\n(%s)", length(h_m$acc))),
#   cat.cex=0.7
# )
# grid.newpage()
# draw.pairwise.venn(
#   area1=length(m_r$acc),
#   area2=length(m_m$acc),
#   cross.area=length(m_overlap_ids),
#   category=c(sprintf("RNA-seq\n(%s)", length(m_r$acc)), 
#              sprintf("compendia\n(%s)", length(m_m$acc))),
#   cat.cex=0.7
# )
# 
# 
# # from here out --> microarray + rnaseq SEPARATELY

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


stopifnot(length(unique(comb_metadata$sample_acc))==nrow(comb_metadata))
#write_csv(comb_metadata, "data/01_metadata/combined_human_mouse_meta.csv")

#comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")
# ----- 3. COUNT TABLES ----- #
# TODO: UPDATE THE COUNTS  / overlap

# // TODO: what is our cutoff ?
cutoff <- comb_metadata %>%
  mutate(expr_sex=case_when(
    is.na(expr_sex) ~ "unknown",
    p_male < 0.7 & p_male > 0.3 ~ "mixed",
    TRUE ~ expr_sex
  ))

ggplot(cutoff, 
       aes(x=data_type, fill=expr_sex)) +
  geom_bar()+
  facet_grid(.~organism)


by_sample <- comb_metadata %>% 
  select(-p_male) %>%
  mutate(expr_sex=ifelse(is.na(expr_sex), "unknown", expr_sex)) %>%
  pivot_longer(cols=c(expr_sex, metadata_sex), 
               names_to="label_type",
               values_to="sex_lab") %>%
  # change names so prettier
  mutate(label_type=case_when(label_type=="metadata_sex" ~ "metadata",
                              label_type=="expr_sex" ~ "expression")) %>%
  mutate(data_type=case_when(data_type=="RNA-seq" ~ "rnaseq",
                             TRUE ~ data_type)) %>%
  mutate(label_type=factor(label_type, levels=c("metadata", "expression")))

# check we havent added samples! (multiplied by 2 b/c separated into expr/metadata)
stopifnot(nrow(by_sample)/2==length(unique(comb_metadata$sample_acc)))
stopifnot(nrow(by_sample)/2==nrow(comb_metadata))
stopifnot(nrow(by_sample)/2==length(unique(by_sample$sample_acc)))
         
ggplot(by_sample %>%
         mutate(sex_lab = ifelse(sex_lab=="mixed", "mixed/pooled", sex_lab)) %>%
         mutate(sex_lab=factor(sex_lab, 
                               levels=c("female", "mixed/pooled", "male", "unknown"))) %>%
         rename(`sample sex`=sex_lab), 
       aes(x=data_type, fill=`sample sex`))+
  geom_bar()+
  facet_grid(label_type ~organism)+  
  xlab("")+
  ylab("Number of samples")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(labels = comma)

ggsave("figures/paper_figs/fig1_sample.png")


# by-study counts - NOTE: slow
by_study <- by_sample %>%
  separate_rows(study_acc, sep=";") %>%
  group_by(study_acc, organism, data_type, label_type) %>%
  summarize(num_tot=n(),
         num_f=sum(sex_lab=="female", na.rm=TRUE),
         num_m=sum(sex_lab=="male", na.rm=TRUE))

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

ggplot(by_study2 %>% 
         rename(`study sex`=study_sex), 
       aes(x=data_type, fill=`study sex`)) +
  geom_bar()+
  facet_grid(label_type~organism)+
  xlab("")+
  ylab("Number of studies")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("figures/paper_figs/fig1_study.png")



# --- create Supplementary Table 1 --- #
# missingness - these tables
#  | organism | dataset | platform | sex | num_studies
#  | organism | dataset | platform | sex | num_samples
counts_sample <- by_sample %>%
  group_by(organism, data_type, label_type, sex_lab) %>%
  count() %>%
  ungroup() 

counts_sample_table <- counts_sample %>% 
  pivot_wider(names_from=sex_lab, values_from=n, values_fill = list(n=0)) %>%
  mutate(num_samples=(female+male+mixed+unknown)) %>%
  mutate(frac_missing=unknown/num_samples,
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


# are the differences significant?
counts_sample_table %>% 
  filter(label_type=="metadata") %>% 
  select(organism, data_type, unknown, num_samples) %>%
  mutate(num_present=num_samples-unknown) %>%
  select(-num_samples) 

# // TODO: stop hard-coding
mat1 <- matrix(c(233526,96982,195258, 34531), byrow=TRUE, nrow=2)
mat2 <- matrix(c(83926 ,39353, 298072,61070), byrow=TRUE, nrow=2)
rownames(mat1) <- c("microarray", "rnaseq")
colnames(mat1) <- c("missing", "present")
human_sample_test <- chisq.test(mat1)

human_sample_test$expected
human_sample_test$observed

# --- create Supplementary Table N --- #
# overlap


# <-- where did I look into this ("05_analysis/01_compendia_overlap.R")






