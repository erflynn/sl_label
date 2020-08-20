# 03_cell_breakdown.R
# E Flynn
# 7/27/2020
#
# Plot the cell line breakdown and sample switching.
#
# Tables created (3)
# Figures created (3)
#
# TODOs:
# - add a figure illustrating sex breakdown separated by tissue
# - cell line cluster densities
# - check on mouse results - are they ok? FIX 
# - extract addtl info on # reports from cellosaurus...
# - divide this file into multiple!


require('tidyverse')

#  --- set up color scheme --- #
library(RColorBrewer)
my.l <- brewer.pal(n = 8, name = 'Set2')
blues <- brewer.pal(9,name="Blues")
oranges <- brewer.pal(9,name="Oranges")
my.cols3 <- c(my.l[3],my.l[2], my.l[8])
my.cols4 <- c(my.l[3],my.l[4],my.l[2],  my.l[8])
my.cols5 <- c(my.l[3], my.l[2], my.l[4], my.l[5], my.l[8])
my.cols6 <- c(my.l[3], blues[4],my.l[4],oranges[4], my.l[2], my.l[8])
# ------- #

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv",col_types="cccccccdld")


# --- CELL LINE LABELS --- #
compendia_cl <- read_csv("data/02_labeled_data/human_compendia_sample_cl.csv")
rnaseq_cl <- read_csv("data/02_labeled_data/human_rnaseq_sample_cl.csv")
# // TODO: aren't there more sample CL files??

m_compendia_cl <- read_csv("data/02_labeled_data/mouse_compendia_sample_cl.csv")
m_rnaseq_cl <- read_csv("data/02_labeled_data/mouse_rnaseq_sample_cl.csv")

cl_lab_unique <- bind_rows(compendia_cl, 
                           rnaseq_cl, 
                           m_rnaseq_cl, 
                           m_compendia_cl) %>% unique() # 63843

# ---- SAMPLE SOURCE TYPE ---- #
options(stringsAsFactors=FALSE)
sample_metadata <- read.csv("data/01_metadata/human_metadata.csv")
rnaseq_sample_metadata <- read.csv("data/01_metadata/human_rnaseq_sample_metadata.csv")
m_sample_metadata <- read.csv("data/01_metadata/mouse_metadata.csv")
m_rnaseq_sample_metadata <- read.csv("data/01_metadata/mouse_rnaseq_sample_metadata.csv")


meta_unique <- sample_metadata %>% 
  select(acc, cl_line, part, title) %>% 
  bind_rows(rnaseq_sample_metadata %>% 
              select(acc, cl_line, part, title)) %>%
  bind_rows(m_sample_metadata %>% 
              select(acc, cl_line, part, title)) %>% 
  bind_rows(m_rnaseq_sample_metadata %>% 
              select(acc, cl_line, part, title)) %>%
  unique() # 560297

mu_s <- meta_unique %>% 
  left_join(cl_lab_unique, by=c("acc"="gsm")) %>%
  group_by(acc) %>%
  mutate(str=paste(c(cl_line, part, title), collapse=" ")) %>%
  ungroup() %>%
  mutate(str=tolower(str))
mu_s2 <- mu_s %>% mutate(source_type=case_when(
  str_detect(str, "xenograft") ~ "xenograft",
  str_detect(str, "stem cell|ipsc") ~ "stem_cell",
  !is.na(accession) & !str_detect(str, "primary") ~ "named_cl",
  !is.na(accession) & str_detect(str, "primary") ~ "other",
  str_detect(str, "tumor|cancer|carcinoma|melanoma|malignant") ~ "cancer",
  str_detect(str, "primary") ~ "primary_cells",
  str_detect(str, "cell line") ~ "unnamed_cl",
  str_detect(str, "cell|culture|passage")  ~ "other",
  is.na(cl_line) | cl_line %in% c("", " ", "--") ~ "tissue",
  
  # // TODO - fix this for mouse -- much of this is the type "c57bl/6, b6, etc"
  TRUE ~ "other"
))

mu_s2 %>% filter(source_type=="other") %>% sample_n(10) %>%
  select(cl_line)

table(mu_s2$source_type)

# Table (1): sample source type
stopifnot(nrow(comb_metadata)==nrow(mu_s2))
mu_s2 %>% write_csv("data/sample_source_type.csv")

# Figure (A): density
comb_metadata_w_src <- comb_metadata %>% 
  inner_join(mu_s2 %>% 
               select(acc, source_type), by=c("sample_acc"="acc"))

tissue_data <- comb_metadata_w_src %>%
  filter(source_type == "tissue" & !is.na(p_male))

stopifnot(nrow(comb_metadata_w_src)==nrow(comb_metadata))

ggplot(comb_metadata_w_src %>% 
         filter(source_type %in% c("named_cl",
                                   "unnamed_cl","tissue","primary_cells",
                                   "stem_cell", "cancer")) %>%
         mutate(source_type=fct_recode(source_type, "cell line"="named_cl", 
                                       "primary cells"="primary_cells",
                                       "stem cells"="stem_cell",
                                       "cancer cells"="cancer")) %>%
         mutate(source_type=fct_collapse(source_type, "cell line"=c("cell line", "unnamed_cl"))) %>%
         rename("sample source"=source_type) %>%
         unite(col="data_src", c("organism", "data_type"), sep=" - "), 
       aes(x=p_male, col=`sample source`)) +
  geom_density()+
  theme_bw()+
  facet_grid(metadata_sex~data_src, scales="free")+
  xlab("P(male)")
ggsave("figures/paper_figs/cl_tiss_density.png")

counts_by_sample_src <- comb_metadata_w_src %>% 
  group_by(organism, data_type, source_type) %>%
  count() %>%
  pivot_wider(names_from=source_type, values_from=n) %>%
  select(organism, data_type, tissue, named_cl, unnamed_cl, primary_cells, stem_cell, cancer, xenograft, tissue) %>%
  rename("cell line (named)"=named_cl,
         "cell line (unnamed)"=unnamed_cl,
         "cancer cells"=cancer)

counts_by_sample_src %>% write_csv("tables/counts_by_sample_src.csv")



# mu_s2 %>% 
#   left_join(comb_metadata_w_src %>% 
#               select(sample_acc, organism, 
#                      data_type, metadata_sex, p_male), by=c("acc"="sample_acc")) %>%
#   filter(source_type=="primary_cells", 
#          organism=="mouse", 
#          p_male > 0.5, p_male < 0.7, 
#          data_type=="rnaseq") %>% 
#   select(acc, part, data_type, metadata_sex, p_male) %>%
#   sample_n(50)


mu_s3 <- mu_s2 %>% 
  select(acc, source_type, cl_line, accession) %>%
  rename(sample_acc=acc, cl_acc=accession, cl_name=cl_line) %>%
  left_join(comb_metadata, by=c("sample_acc")) %>%
  select(sample_acc, organism, data_type, everything()) 

### -- cell line vs tissue sex breakdown -- ##
mu_s2 <- read_csv("data/sample_source_type.csv")

sex_lab_w_source <-mu_s2 %>% 
  select(acc, source_type) %>% 
  rename(sample_acc=acc) %>%
  left_join(comb_metadata) %>%
  select(-platform)


frac_dat <- sex_lab_w_source %>% 
  select(source_type, data_type, organism, expr_sex) %>%
  filter(!is.na(expr_sex)) %>%
  mutate(source_type2=case_when(
    source_type=="tissue" ~ "tissue",
    source_type %in% c("named_cl", "unnamed_cl") ~ "cell line",
    source_type=="primary_cells" ~ "primary cells",
    TRUE ~ "other"
  )) %>%
  bind_rows(sex_lab_w_source %>% 
              select(source_type, data_type, organism, expr_sex) %>%
              filter(!is.na(expr_sex)) %>%
              mutate(source_type2="aggregate")) %>%
  group_by(data_type, source_type2, organism) %>% 
  mutate(total_f=sum(expr_sex=="female"),
         total_m=sum(expr_sex=="male")) %>%
  ungroup() %>%
  group_by(source_type2, organism, data_type) %>%
  mutate(num_f=sum(expr_sex=="female"),
         num_m=sum(expr_sex=="male")) %>%
  mutate(frac_f=num_f/(total_f+total_m),
         frac_m=num_m/(total_f+total_m)) %>%
  select(-source_type, -expr_sex) %>%
  unique() %>% 
  ungroup()

frac_dat2 <- frac_dat %>% 
  select(data_type, organism,  source_type2, total_f, total_m, frac_f, frac_m) %>%
  pivot_longer(cols=c("frac_f", "frac_m"), 
               names_to="sex", 
               values_to="fraction") %>%
  pivot_longer(cols=c("total_f", "total_m"), 
               names_to="sex2", 
               values_to="total") %>%
  mutate(sex=str_replace(sex, "frac_", ""),
         sex2=str_replace(sex2, "total_", "")) %>%
  filter(sex==sex2) %>%
  select(-sex2) %>%
  pivot_longer(cols=c("total", "fraction"), 
               names_to="stat",
               values_to="value")

require('scales')
ggplot(frac_dat2 %>% 
         mutate(organism=sprintf("%s - %s", organism, data_type)) %>%
         mutate(source_type2=factor(source_type2, 
                                    levels=c("aggregate", "cell line",
                                             "tissue", "other"))) %>%
         mutate(sex=ifelse(sex=="m", "male", "female")), 
       aes(x=source_type2, y=value, fill=sex))+
  geom_bar(stat="identity")+facet_grid(stat~organism, scales="free")+theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.title = element_blank()) +
  scale_fill_manual(values=c(my.cols4[1], my.cols4[3]))+
  scale_y_continuous(labels = comma)+
  xlab("")+
  ylab("")
ggsave("figures/paper_figs/sex_breakdown_sample_type.png")

# breakdown calculations for stats
frac_dat2 %>% 
  filter(sex=="f" & source_type2 %in% c("tissue", "aggregate") & stat=="fraction") %>%
  arrange(organism, data_type, source_type2)

frac_dat2 %>% 
  filter(source_type2 %in% c("cell line", "aggregate"), stat=="total") %>%
  arrange(organism, data_type, source_type2) %>% 
  group_by(organism, data_type, source_type2) %>% 
  summarize(sum=sum(value))



# distributions are different
ggplot(sex_lab_w_source %>%
         filter(!is.na(expr_sex)),
        aes(x=p_male, group=source_type, col=source_type))+
   geom_density()+
  facet_wrap(.~organism) # // TODO: mouse "other" looks weird


# --- CELL LINE SEX LABELS --- #

cell_df <- read.csv("data/00_db_data/cellosaurus_df_v2.txt") 
# // TODO v3 breaks alleles, use allele_freq3
cell_sex <- cell_df %>% 
  select(primary_accession, cl, sex, alleles) %>% 
  mutate(sex=tolower(sex))  %>%
  mutate(sex=case_when(
    sex %in% c("", "sex unspecified") ~ "unknown" ,
    sex %in% c("mixed sex","sex ambiguous") ~ "mixed",
    TRUE ~ sex
  )) %>%
  mutate(allele_sex=case_when(
    alleles %in% c("", "Not_detected") ~ "unknown",
    alleles %in% c("X,Y|X", "X|X,Y") ~ "both",
    str_detect(alleles, "Y") ~ "male",
    str_detect(alleles, "X") ~ "female",
    TRUE ~ "unknown"
  )) %>%
  rename(annot_sex=sex) %>%
  mutate(primary_accession=tolower(primary_accession))

# Table (1) - cell line sex labels
# //TODO: move this
cell_sex %>% write_csv("data/ref_cell_line_sex_labels.csv")

cell_sex2 <- cell_sex %>%
  rename(cl_acc=primary_accession, 
         cl_name=cl,
         cl_annot_sex=annot_sex,
         cl_allele_sex=allele_sex) %>%
  select(-alleles) %>%
  unique()

stopifnot(length(unique(cell_sex2$cl_acc))==nrow(cell_sex2))

# separate out multi-mapping
cl_acc_maps <- mu_s3 %>% 
  filter(!is.na(cl_acc)) %>% 
  select(cl_acc) %>% 
  unique() 
multi_cl_acc <- cl_acc_maps %>% 
  filter(str_detect(cl_acc, ";")) %>%
  mutate(cl_acc_split=cl_acc) %>%
  mutate(id=1:n()) %>%
  separate_rows(cl_acc_split, sep=";") %>%
  unique() %>%
  arrange(id, cl_acc_split)

multi_cl_acc_w_lab <- multi_cl_acc %>% 
  left_join(cell_sex2, by=c(cl_acc_split="cl_acc")) %>%
  #mutate(cl_annot_sex=ifelse(cl_annot_sex=="unknown", NA, cl_annot_sex),
  #       cl_allele_sex=ifelse(cl_allele_sex=="unknown", NA, cl_allele_sex)) %>%
  group_by(id) %>%
  arrange(cl_acc_split) %>%
  mutate(cl_acc_split=paste(cl_acc_split, collapse=";"),
         cl_name=paste(cl_name, collapse=";"),
         cl_annot_sex=paste(sort(unique(cl_annot_sex[!is.na(cl_annot_sex)])), collapse=";"),
         cl_allele_sex=paste(sort(unique(cl_allele_sex[!is.na(cl_allele_sex)])), collapse=";")) %>%
  unique() %>%
  ungroup() %>%
  select(-id) 

single_map <- cl_acc_maps %>% filter(!str_detect(cl_acc, ";")) %>%
  left_join(cell_sex2, by="cl_acc")  %>%
  mutate(cl_acc_split=cl_acc) %>%
  select(colnames(multi_cl_acc_w_lab))

cl_lab_comb <- multi_cl_acc_w_lab %>% bind_rows(single_map) 
stopifnot(nrow(cl_lab_comb)==length(unique(cl_lab_comb$cl_acc)))

samples_cl_sl <- mu_s3 %>% left_join(cl_lab_comb %>% select(-cl_name), by="cl_acc") %>%
  mutate(cl_acc_split=cl_acc) %>% # split column has things properly formatted
  select(-cl_acc_split)
stopifnot(length(unique(samples_cl_sl$sample_acc))==nrow(samples_cl_sl))

# Table (2) - samples mapped to cell lines with sex labels
samples_cl_sl %>% write_csv("data/cl_sex_mapped.csv")
samples_cl_sl <- read_csv("data/cl_sex_mapped.csv")

cl_expr <- samples_cl_sl %>% 
  filter(cl_annot_sex %in% c("female", "male") & !is.na(expr_sex)) %>%
  mutate(tot=n()) %>%
  group_by(organism,cl_annot_sex, expr_sex) %>%
  mutate(n=n()) %>%
  select(cl_annot_sex, expr_sex, n, tot) %>%
  unique() %>%
  mutate(frac=n/tot) %>%
  ungroup()
cl_expr # <-- this is the confusion matrix for percents that switch!


# -- Figure (1) -- #
# alluvial diagram for cell line sex switching 
require("ggalluvial")
my_dat3 <- samples_cl_sl %>% 
  filter(source_type=="named_cl") %>%
  rename(allele_sex=cl_allele_sex, annot_sex=cl_annot_sex) %>%
  mutate(
  allele_sex=case_when(
    (allele_sex %in% c("female")) ~ "female",
    (allele_sex %in% c("male")) ~ "male",
    allele_sex=="both" ~ "both",
    allele_sex=="unknown" ~ "unknown",
    TRUE ~ "multi-map"
    
  )) %>%
  mutate(
    annot_sex=case_when(
      (annot_sex %in% c("female")) ~ "female",
      (annot_sex %in% c("male")) ~ "male",
      annot_sex=="both" ~ "both",
      annot_sex=="unknown" ~ "unknown",
      TRUE ~ "multi-map"
    )) %>%
  filter(!is.na(expr_sex) & !is.na(allele_sex) & !is.na(annot_sex)) %>%
  select(-cl_acc,-metadata_sex, -cl_name, -organism, -platform,
         -source_type, -study_acc, -p_male, -data_type) 
  


cl_flow_freq_counts <- my_dat3 %>%
  ungroup() %>%
  filter(!is.na(expr_sex) & expr_sex!="") %>%
  mutate(allele_sex=ifelse(is.na(allele_sex) | allele_sex=="", "unknown", allele_sex),
         annot_sex=ifelse(is.na(annot_sex) | annot_sex=="", "unknown", annot_sex )) %>%
  rename(donor=annot_sex, expression=expr_sex, recorded=allele_sex) %>%
  select(-sample_acc) %>%
  group_by(donor, expression, recorded) %>% 
  mutate(Freq=n()) %>% 
  unique() %>%
  ungroup() %>%
  mutate(row_id=1:n()) %>%
  gather(key="labeling_method", value="sex", -Freq, -row_id) %>%
  mutate(row_id=as.factor(row_id), 
         labeling_method=factor(labeling_method, levels=c("donor", "recorded", "expression")),
         sex=as.factor(sex)) %>%
  unique() %>%
  mutate(sex=factor(sex, 
                    levels=c("female","male", 
                             "both", "multi-map", "unknown")))

ggplot(cl_flow_freq_counts,
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
  theme(legend.position = "none")+
  scale_fill_manual(values=my.cols5)

ggsave("figures/paper_figs/supp_fig_alluv_donor_recorded.png")

# simpler alluv fig
cl_flow_freq_counts2 <- my_dat3 %>%
  ungroup() %>%
  select(-allele_sex) %>%
  filter(!is.na(expr_sex) & expr_sex!="" ) %>%
  mutate(annot_sex=ifelse(is.na(annot_sex) | annot_sex=="", "unknown", annot_sex)) %>%
  rename(donor=annot_sex, expression=expr_sex) %>%
  select(-sample_acc) %>%
  group_by(donor, expression) %>% 
  mutate(Freq=n()) %>% 
  unique() %>%
  ungroup() %>%
  mutate(row_id=1:n()) %>%
  gather(key="labeling_method", value="sex", -Freq, -row_id) %>%
  mutate(row_id=as.factor(row_id), 
         labeling_method=factor(labeling_method, levels=c("donor", "expression")),
         sex=as.factor(sex)) %>%
  unique() %>%
  mutate(sex=factor(sex, 
                    levels=c("female","male", "multi-map", "unknown")))

ggplot(cl_flow_freq_counts2,
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
  theme(legend.position = "none")+
  scale_fill_manual(values=c(my.cols5[1], my.cols5[2],my.cols5[4], my.cols5[5]))
ggsave("figures/paper_figs/fig_alluv_cl_simple.png")

# --- rearrange the sample level labels to cell line level labels -- #
# add columns for
#  - number of female or male samples per cl line
#  - avg/se for sex labeling score
#  - condensed expression sex labels (e.g. "female" or "female;male")
#  - condensed annotation labels
#  - condensed reported sex
cl_level_lab <- samples_cl_sl %>% 
  filter(!is.na(cl_acc)) %>%
  select(-sample_acc, -organism, 
         -data_type, -source_type, 
         -metadata_sex, -study_acc, 
         -platform) %>%
  group_by(cl_acc) %>%
  mutate(num_f=sum(expr_sex=="female", na.rm=TRUE), 
         num_m=sum(expr_sex=="male", na.rm=TRUE),
         num_unknown=sum(is.na(expr_sex)),
         avg_p_male=mean(p_male, na.rm=TRUE), 
         se_p_male=sd(p_male, na.rm=TRUE),
         expr_sex=paste(sort(unique(expr_sex[!is.na(expr_sex)])), collapse=";"),
         cl_annot_sex=paste(sort(unique(cl_annot_sex[cl_annot_sex!="unknown"])), collapse=";"),
         cl_allele_sex=paste(sort(unique(cl_allele_sex[cl_allele_sex!="unknown"])), collapse=";")) %>%
  select(-p_male, -cl_name) %>%
  unique() %>%
  mutate(frac_f=num_f/(num_m+num_f),
         frac_m=num_m/(num_m+num_f),
         frac_unknown=num_unknown/(num_m+num_f)) %>%
  ungroup() %>%
  unique()
cl_level_lab %>% 
  filter(!str_detect(cl_acc, ";")) %>% 
  arrange(cl_annot_sex, frac_f) 

# remove the multi-mapping to take a quick look
single_map <- cl_level_lab %>% 
  filter(!str_detect(cl_acc, ";"))
table(single_map$cl_annot_sex, single_map$expr_sex)

stopifnot(length(unique(cl_level_lab$cl_acc))==nrow(cl_level_lab))


# add additional columns to the cell line level data
#  - previously documented changes (`annot_all`)
#  - number of samples
#  - 95% ci 
df <- cl_level_lab %>% 
  filter(!str_detect(cl_acc, ";")) %>%
  mutate(cl_annot_sex=ifelse(cl_annot_sex=="", 
                             "unknown",cl_annot_sex),
         cl_allele_sex=ifelse(cl_allele_sex=="", 
                              "unknown",cl_allele_sex)) %>%
  filter(cl_annot_sex!="unknown" & expr_sex!="") %>%
  group_by(cl_acc) %>%
  mutate(annot_all=paste(cl_annot_sex,cl_allele_sex, sep=">")) %>%
  ungroup() %>%
  mutate(num_s=num_f+num_m) %>%
   select( cl_acc, 
          cl_annot_sex, cl_allele_sex, expr_sex, annot_all,
          num_s, frac_f, frac_m, avg_p_male, se_p_male) %>%
  unique() %>%
  mutate(
    ci_l=(avg_p_male-1.96*se_p_male),
    ci_u=(avg_p_male+1.96*se_p_male)
  )

stopifnot(length(unique(df$cl_acc  ))==nrow(df))
table(df$annot_all, df$expr_sex)

  

# high confidence switching
# We define high confidence switching as a 95% CI for 
# cell line scores > 0.6 or < 0.4.
df2 <- df %>%
  group_by(annot_all) %>%
  mutate(n=n()) %>%
  mutate(num_s=ifelse(num_s > 100, 100, num_s)) %>%
  rename(`number of samples`=num_s,
         `annotated sex`=cl_annot_sex) %>%
  ungroup() %>%
  mutate(annot_all=sprintf("%s\n(n=%s)", annot_all, n))


# Figure (2) #
# ---- SUPPLEMENTARY FIGURE A ---- #
# how does avg expression score correspond to to donor sex
# and reported sex?

#  male>male & p_male < 0.5
#  male>unknown & p_male < 0.5
#  female>unknown & p_male > 0.5
#  female>female & p_male > 0.5
ggplot(df2 %>%
         mutate(
           `ychr loss`=factor(case_when(
           (cl_allele_sex=="male" & `annotated sex`=="male" & avg_p_male < 0.5) ~ "novel",
           (cl_allele_sex=="unknown" & `annotated sex`=="male" & avg_p_male < 0.5) ~ "novel",
           (cl_allele_sex=="female" & `annotated sex`=="male" & avg_p_male < 0.5) ~ "documented",
           (cl_allele_sex=="both" & `annotated sex`=="male" & avg_p_male < 0.5 )~ "documented",
           TRUE ~ "N/A"
         ), levels=c("documented", "novel", "N/A"))),
       aes(y=avg_p_male, x=1))+
  geom_boxplot(outlier.alpha=0)+
  geom_point(aes(size=`number of samples`, col=`ychr loss`),
             position=position_jitter(0.5), 
             alpha=0.1)+
  scale_size(breaks=c(1,5,10,25, 50, 100))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylim(c(0,1))+
  ylab("Average P(male)")+
  xlab("donor sex > reported sex")+
  geom_hline(yintercept=0.5, lty=2, col="gray")+
  facet_grid(.~annot_all)+
  scale_color_manual(values=c("blue", "purple", "black"))

ggsave("figures/paper_figs/supp_cl_reported_donor.png")
# // TODO: 
#  - deal w tiny categories


# we found novel data
#  male>male & p_male < 0.5
#  male>unknown & p_male < 0.5
#  female>unknown & p_male > 0.5
#  female>female & p_male > 0.5

# add a cell line switching category
#  - hc_switch: switching with 95% CI > 0.6 or < 0.4
#  - hc_no_switch: no switching with 95% CI > 0.6 or < 0.4
#  - some_switch: more than one sample of the alternate sex
df_switch <- df %>%
  mutate(switching_category=case_when(
    num_s <= 3 ~ "insufficient samples",
    (cl_annot_sex=="female" & ci_l > 0.6) | 
      (cl_annot_sex=="male" & ci_u < 0.4) ~ "hc_switch",
    (cl_annot_sex=="female" & ci_u < 0.4) | 
     (cl_annot_sex=="male" & ci_l > 0.6) ~ "hc_no_switch",
     ((cl_annot_sex=="female" & frac_m*num_s > 1) |
        (cl_annot_sex=="male" & frac_f*num_s > 1))  ~ "some_switch",
    TRUE ~ "no_switch"
  )) %>%
  mutate(doc_switch=case_when(
    annot_all %in% c("female>both" ,"female>male", "male>female", "male>both") ~ "doc_switch",
    annot_all %in% c("female>female", "male>male") ~ "doc_no_switch",
    TRUE ~ "unknown"
  ))
df_switch2 <- df_switch %>%
  filter(switching_category!="insufficient samples")

table(df_switch2$switching_category, df_switch2$doc_switch)

table(df_switch2$switching_category, df_switch2$cl_annot_sex)
table(df_switch2$switching_category, df_switch2$cl_annot_sex)/nrow(df_switch2)



# --- use a diptest to separate the data into uni-modal and multi-modal --- #
# // TODO: really should be doing clustering instead of dip-test
require('diptest')

# add a column with the dip-test p-value
dps <-samples_cl_sl %>% 
  filter(!is.na(cl_acc)) %>% 
  group_by(cl_acc) %>%
  mutate(dip_p=dip.test(p_male)$p)

# look specifically at the some-switch data with 40-60% of of each sex
# add a column for dip_pass (TRUE == unimodal; null hypothesis is unimodal)
dps_some_swtch <- dps %>% 
  filter(!is.na(cl_acc)) %>% 
  select(cl_name, p_male, dip_p) %>% 
  right_join(df_switch2 %>% 
               filter(frac_m > 0.40 & frac_m < 0.60),
             by="cl_acc") %>%
  filter(switching_category=="some_switch") %>% 
  mutate(dip_pass=ifelse((dip_p > 0.05), "unimodal", "multimodal"))
# TODO -- do we need to multiple hypothesis correct here?

# breakdown uni- vs multi-modal
table(dps_some_swtch %>% select(cl_acc, dip_pass) %>% 
        unique() %>% pull(dip_pass))


# -- plot what happens with some switching data -- #
# // TODO normalize density
dps_some_swtch %>%
  rename(`Cell line`=cl_acc) %>%
  ggplot(aes(x=p_male, col=`Cell line`))+geom_density()+
  facet_grid(dip_pass~., scales="free")+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("P(male)")

dip_test_col <- dps %>% select(cl_acc, dip_p) %>% unique() %>% ungroup()
stopifnot(length(unique(dip_test_col$cl_acc))==nrow(dip_test_col))

# create and save a table
# - add the dip-test info
# - add cell name 
# - rename columns
df_switch3 <- df_switch2 %>% 
  left_join(cell_sex2 %>% select(cl_acc, cl_name)) %>%
  left_join(dip_test_col) %>%
  select(cl_acc, cl_name, switching_category, cl_annot_sex, cl_allele_sex, 
         avg_p_male, se_p_male,  num_s, frac_f, frac_m, dip_p) %>% 
  mutate(across(contains("frac|p_male|ci|dip_p"), signif, 3)) %>%
  rename(`Accession`=cl_acc, `Name`=cl_name, `Category`=switching_category,
         `Annotated sex`=cl_annot_sex, `Reported sex (amel.)`=cl_allele_sex,
         `Average P(male)`=avg_p_male, `SE P(male)`=se_p_male,
         `Number of samples`=num_s, `Expr. fraction female`=frac_f, 
         `Expr. fraction male`=frac_m,`Dip-test Pval`=dip_p) %>%
  arrange(Category,`Average P(male)`)
  
stopifnot(nrow(df_switch3)==nrow(df_switch2))
# Table (3) - cell lines + switching category + additional info
df_switch3 %>%
  write_csv("tables/supp_cl_sample_switch.csv")


# Figure (3) #
# ---- HC switch figure ---- #
# add a "novel" column to denote if it's new
hc_switch <- df_switch %>% 
  filter(switching_category=="hc_switch") %>%
  select(-expr_sex, -frac_f, -frac_m, -se_p_male) %>%
  mutate(novel=case_when(cl_allele_sex == "unknown" |
                   cl_allele_sex==cl_annot_sex ~ "novel",
         TRUE ~ "documented")) %>%
  arrange(avg_p_male) %>% 
  left_join(cell_sex2 %>% select(cl_acc, cl_name)) %>%
  mutate(cl_name=factor(cl_name, levels=cl_name)) 

table(hc_switch$novel, hc_switch$cl_annot_sex)

# // TODO: geom_label_repel?? for the novel ones?
ggplot(hc_switch %>%
         mutate(num_s=ifelse(num_s > 100, 100, num_s)) %>%
         rename(`number of samples`=num_s), 
       aes(y=avg_p_male, x=cl_name, col=novel))+
  scale_fill_manual(values=c("black", "purple"))+
  scale_color_manual(values=c("black", "purple"))+
  geom_point(alpha=0.2, aes(size=`number of samples`))+
  scale_size(breaks=c(1,5,10,25, 50, 100))+
  geom_errorbar(aes(ymin=ci_l, ymax=ci_u), 
                alpha=0.2,width=0.3)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 90, size=7,
                                   vjust = 0.5, hjust=1))+
  ylab("Average P(male)")+
  xlab("Cell line")+
  #geom_text(aes(label=cl_name), size=2, 
  #              position=position_jitter(0.1))+
  geom_hline(yintercept=0.5, lty=2, col="gray")+
  geom_hline(yintercept=1, lty=1, col="gray")+
  geom_hline(yintercept=0, lty=1, col="gray")



# ---- allele frequency in text ---- #
allele_freq <- read_csv("data/00_db_data/cellosaurus_allele_freq.csv")

allele_freq2 <- allele_freq %>% filter(`Not_detected`==0 & Y==0) %>%
  select(-Y, -Not_detected) %>%
  mutate(across(c(`X,Y`,`X`), ~./num_srcs)) %>%
  arrange(`X`)

ggplot(allele_freq2 %>% 
         mutate(cl_acc=factor(cl_acc, levels=unique(allele_freq2$cl_acc))), 
       aes(x=cl_acc, y=`X,Y`, size=num_srcs))+
  geom_point(alpha=0.5)+
  ylab("fraction XY chromosome")+
  xlab("cell line")+
  theme_bw()+
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_blank())

# allele_freq2 %>%
#   mutate(cl_acc=factor(cl_acc, levels=unique(allele_freq2$cl_acc))) %>%
#   ggplot(aes(x=cl_acc, y=frac_src, fill=allele))+geom_bar(stat="identity")

compare_src <- allele_freq2 %>% 
  mutate(cl_acc=tolower(cl_acc)) %>%
  inner_join(cl_level_lab , by=c("cl_acc")) %>%
  select(cl_acc, num_srcs, `X,Y`, `X`, `frac_f`, `frac_m`, `avg_p_male`, `se_p_male`) %>%
  filter(!is.na(frac_m)) %>%
  left_join(dip_test_col) %>%
  mutate(dip_pass=ifelse((dip_p > 0.05), "unimodal", "multimodal")) %>%
  mutate(
    ci_l=(avg_p_male-1.96*se_p_male),
    ci_u=(avg_p_male+1.96*se_p_male)
  )

cor.test(compare_src$`X,Y`, compare_src$avg_p_male) # 0.218


ggplot(compare_src %>%
         rename("Number of STR refs"=num_srcs,
                "SE p-male"=se_p_male), 
       aes(x=`X,Y`, y=avg_p_male, col=`Number of STR refs`, size=`SE p-male`))+
  geom_point(alpha=0.5)+
  ylim(0,1)+
  #geom_errorbar(aes(ymin=ci_l, ymax=ci_u))+
  facet_wrap(.~dip_pass)+
  theme_bw()+
  ylab("Average p-male")+
  xlab("Fraction male STR reported")
ggsave("figures/paper_figs/compare_str_fraction.png")

# // TODO - is it worth adding the reference breakdown?


# -------- LOOK AT STUDIES ------- #
# is it between study heterogeneity or within study heterogeneity?
head(hc_switch)

study_grp <- samples_cl_sl %>%
  semi_join(df_switch, by="cl_acc") %>%
  separate_rows(study_acc, sep=";") %>%
  group_by(cl_acc) %>%
  group_by(cl_acc, study_acc) %>%
  summarize(num_samples=n(),
          avg_p_study=mean(p_male),
          se_p_study=sd(p_male)) %>%
  ungroup() %>%
  group_by(cl_acc) %>%
  mutate(num_studies=n()) %>%
  ungroup()

cls <- study_grp %>% filter(num_studies > 3, num_samples > 3) %>% 
  distinct(cl_acc) %>% pull()

study_grp %>% 
  left_join(df_switch) %>% 
  left_join(dip_test_col) %>%
  mutate(dip_pass=ifelse((dip_p > 0.05), "unimodal", "multimodal")) %>%
  mutate(cl_acc=factor(cl_acc, levels=df_switch$cl_acc[order(df_switch$avg_p_male)])) %>%
  filter(cl_acc %in% cls &
           switching_category != "insufficient samples") %>%
  ggplot(aes(y=cl_acc, x=avg_p_study, 
             size=num_samples))+
  geom_point(alpha=0.5)+
  theme(axis.ticks.y=element_blank(),
        axis.text.y = element_blank())+
  facet_grid(switching_category ~ dip_pass, scales="free_y")

#### - mislabeled f-->m

mislab <-samples_cl_sl %>% filter(cl_annot_sex=="female", p_male > 0.9, cl_name!="huvec") %>% 
  arrange(cl_allele_sex, desc(p_male))%>% 
  select(-source_type, -platform, -metadata_sex, -expr_sex, -cl_annot_sex) # 589
mislab %>% group_by(cl_acc) %>% count() %>% arrange(desc(n))

# cvcl_l909 -- likely wrong
# try to see if these look different

#### empirical distribution clustering 