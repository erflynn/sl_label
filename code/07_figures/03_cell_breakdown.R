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
# - alluvial diagram w/o amelogenin
# - ADD MOUSE
# - add a figure illustrating sex breakdown separated by tissue
# - cell line cluster densities
# - extract addtl info on # reports from cellosaurus...
# - divide this file into multiple!


require('tidyverse')

organism <- "human"
comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")
human_metadata <- comb_metadata %>% filter(organism=="human") %>% select(-platform)

# --- CELL LINE LABELS --- #
compendia_cl <- read_csv("data/02_labeled_data/human_compendia_sample_cl.csv")
rnaseq_cl <- read_csv("data/02_labeled_data/human_rnaseq_sample_cl.csv")
# // TODO: aren't there more sample CL files??

cl_lab_unique <- bind_rows(compendia_cl, rnaseq_cl) %>% unique() # 63843

# ---- SAMPLE SOURCE TYPE ---- #
sample_metadata <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", organism))
rnaseq_sample_metadata <- read.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", organism))
meta_unique <- sample_metadata %>% select(acc, cl_line, part, title) %>% bind_rows(
  rnaseq_sample_metadata %>% select(acc, cl_line, part, title)) %>%
  unique() # 560297

mu_s <- meta_unique %>% left_join(cl_lab_unique, by=c("acc"="gsm")) %>%
  group_by(acc)  %>% mutate(str=paste(c(cl_line, part, title), collapse=" ")) %>%
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
  TRUE ~ "other"
))

table(mu_s2$source_type)

# Table (1): sample source type
mu_s2 %>% write_csv("data/sample_source_type.csv")
stopifnot(nrow(human_metadata)==nrow(mu_s2))

mu_s3 <- mu_s2 %>% 
  select(acc, source_type, cl_line, accession) %>%
  rename(sample_acc=acc, cl_acc=accession, cl_name=cl_line) %>%
  left_join(human_metadata, by=c("sample_acc")) %>%
  select(sample_acc, organism, data_type, everything()) 

### -- cell line vs tissue sex breakdown -- ##
sample_type <- read_csv("data/sample_source_type.csv")
human_s <- comb_metadata %>% filter(organism=="human")
sex_lab_w_source <-sample_type %>% 
  select(acc, source_type) %>% 
  rename(sample_acc=acc) %>%
  left_join(human_s) %>%
  select(-organism,-platform)


frac_dat <- sex_lab_w_source %>% 
  select(source_type, data_type, expr_sex) %>%
  filter(!is.na(expr_sex)) %>%
  mutate(source_type2=case_when(
    source_type=="tissue" ~ "tissue",
    source_type %in% c("named_cl", "unnamed_cl") ~ "cell line",
    source_type=="primary_cells" ~ "primary cells",
    TRUE ~ "other"
  )) %>%
  group_by(data_type) %>% 
  mutate(total_f=sum(expr_sex=="female"),
         total_m=sum(expr_sex=="male")) %>%
  ungroup() %>%
  group_by(source_type2, data_type) %>%
  mutate(num_f=sum(expr_sex=="female"),
         num_m=sum(expr_sex=="male")) %>%
  mutate(frac_f=num_f/(total_f+total_m),
         frac_m=num_m/(total_f+total_m)) %>%
  select(-source_type, -expr_sex) %>%
  unique() %>% 
  ungroup()

counts_w_aggreg <- frac_dat %>% 
  select(-frac_f, -frac_m) %>%
  pivot_longer(cols=c("total_f", "num_f")) %>%
  mutate(dat_type=ifelse(name=="total_f", "aggregate", "separate")) %>%
  rename("num_f"=value) %>%
  select(-name) %>%
  pivot_longer(cols=c("total_m", "num_m")) %>%
  mutate(dat_type2=ifelse(name=="total_m", "aggregate", "separate")) %>%
  rename("num_m"=value) %>%
  select(-name) %>%
  filter(dat_type==dat_type2) %>%
  select(-dat_type2) %>%
  pivot_longer(cols=c("num_f", "num_m"), names_to="expr_sex", values_to="count") %>%
  mutate(comb_col=paste(source_type2, dat_type, sep=","))

counts_w_aggreg2 <- counts_w_aggreg %>%
  filter(str_detect(comb_col, "separate") |comb_col=="tissue,aggregate") %>%
  mutate(comb_col=str_replace_all(comb_col, ",separate", "")) %>%
  mutate(comb_col=ifelse(comb_col=="tissue,aggregate", "aggregate", comb_col)) %>%
  select(-source_type2, -dat_type) %>%
  mutate(source_type=factor(comb_col, levels=c("aggregate", "cell line", "tissue", "primary cells", "other"))) 


# # Figure (0) - breakdown
# ggplot(sex_lab_w_source %>% 
#          filter(source_type %in% 
#                   c("tissue", "named_cl") &
#                   !is.na(expr_sex)),
#        aes(x=source_type))+
#   geom_bar(aes(fill=expr_sex))+
#   facet_grid(data_type~.)
# 
# # distributions
# ggplot(sex_lab_w_source %>%
#          filter(source_type %in% c("tissue", "named_cl", "unnamed_cl", "primary_cells")), 
#        aes(x=p_male, group=source_type, col=source_type))+
#   geom_density()


# --- CELL LINE SEX LABELS --- #
cell_df <- read.csv("data/00_db_data/cellosaurus_df_v2.txt")
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

cl_expr <- samples_cl_sl %>% 
  filter(cl_annot_sex %in% c("female", "male") & !is.na(expr_sex)) %>%
  mutate(tot=n()) %>%
  group_by(cl_annot_sex, expr_sex) %>%
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
  select(-cl_acc,-metadata_sex, -cl_name, -organism, 
         -source_type, -study_acc, -p_male, -data_type) 
  


flow_freq_counts <- my_dat3 %>%
  ungroup() %>%
  filter(!is.na(expr_sex) & expr_sex!="" ) %>%
  mutate(allele_sex=ifelse(is.na(allele_sex) | allele_sex=="", "unknown", allele_sex),
         annot_sex=ifelse(is.na(annot_sex) | annot_sex=="", "unknown", annot_sex )) %>%
  rename(metadata=annot_sex, expression=expr_sex, amelogenin=allele_sex) %>%
  select(-sample_acc) %>%
  group_by(metadata, expression, amelogenin) %>% 
  mutate(Freq=n()) %>% 
  unique() %>%
  ungroup() %>%
  mutate(row_id=1:n()) %>%
  gather(key="labeling_method", value="sex", -Freq, -row_id) %>%
  mutate(row_id=as.factor(row_id), 
         labeling_method=factor(labeling_method, levels=c("metadata", "amelogenin", "expression")),
         sex=as.factor(sex)) %>%
  unique() %>%
  mutate(sex=factor(sex, levels=c("female","male", "multi-map", "unknown", "both")))

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
  theme(legend.position = "none") 

# ggsave()


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
         -metadata_sex, -study_acc) %>%
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
ggplot(df2,
       aes(y=avg_p_male, x=1))+
  geom_boxplot(outlier.alpha=0)+
  geom_point(aes(size=`number of samples`, col=`annotated sex`), 
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
  facet_grid(.~annot_all)
# number of cell lines in each
# // TODO: 
#  - highlight validation vs novelty?
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
  mutate(dip_pass=ifelse((dip_p > 0.05), "unimodal", "mulitmodal"))
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
  mutate(novel=(cl_allele_sex == "unknown" |
                   cl_allele_sex==cl_annot_sex)) %>%
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

