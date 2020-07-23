# Plot the cell line breakdown and sample switching

# table:
#  sample | metadata_sex | expr_sex | cell_line(y/n) | cl_acc | cl_name | cl_annot_sex | cl_allele_sex |

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

# SAVE THIS FILE!!!
mu_s2 %>% write_csv("data/sample_source_type.csv")
stopifnot(nrow(human_metadata)==nrow(mu_s2))

mu_s3 <- mu_s2 %>% 
  select(acc, source_type, cl_line, accession) %>%
  rename(sample_acc=acc, cl_acc=accession, cl_name=cl_line) %>%
  left_join(human_metadata, by=c("sample_acc")) %>%
  select(sample_acc, organism, data_type, everything()) 




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
cl_expr # <-- this is the confusion matrix for percents


# -- PLOT -- #

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

###
#SAVE
# cell line level
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

single_map <- cl_level_lab %>% 
  filter(!str_detect(cl_acc, ";"))

table(single_map$cl_annot_sex, single_map$expr_sex)

stopifnot(length(unique(cl_level_lab$cl_acc))==nrow(cl_level_lab))



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
# cell line > 0.6 or < 0.4.

df2 <- df %>%
  group_by(annot_all) %>%
  mutate(n=n()) %>%
  mutate(num_s=ifelse(num_s > 100, 100, num_s)) %>%
  rename(`number of samples`=num_s,
         `annotated sex`=cl_annot_sex) %>%
  ungroup() %>%
  mutate(annot_all=sprintf("%s\n(n=%s)", annot_all, n))

ggplot(df2,
       aes(y=avg_p_male, x=1))+
  geom_boxplot(outlier.alpha=0)+
  geom_point(aes(size=`number of samples`, col=`annotated sex`), 
             position=position_jitter(0.5), 
             alpha=0.1)+
  scale_size(breaks=c(1,5,10,25, 50, 100))+
  # geom_errorbar(aes(ymin=avg_p_male-se_p_male, 
  #                   ymax=avg_p_male+se_p_male,
  #                   col=cl_annot_sex), 
  #               position=position_jitter(0.5), 
  #               width=0.01,
  #               alpha=0.1)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylim(c(0,1))+
  ylab("Average P(male)")+
  xlab("annotated sex > amelogenin sex")+
  geom_hline(yintercept=0.5, lty=2, col="gray")+
  facet_grid(.~annot_all)
# number of cell lines in each



# interesting: 
#  male>male & p_male < 0.5
#  male>unknown & p_male < 0.5
#
#  female>unknown & p_male > 0.5
#  female>female & p_male > 0.5

# add a SWITCH CATEGORY
df_switch = df %>%
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
# look at the distribtuion of "some_switch"

samples_cl_sl %>% filter(!is.na(cl_acc)) %>% select(cl_acc, p_male) %>% 
  right_join(df_switch2 %>% select(cl_acc, num_s, frac_f, cl_annot_sex, 
                                   switching_category, doc_switch), by="cl_acc") %>%
  filter(switching_category=="some_switch") %>%
  ggplot(aes(x=p_male))+geom_histogram()+facet_grid(cl_annot_sex~.)

df_switch2 %>% filter(switching_category=="some_switch") %>%
  ggplot(aes(x=frac_f))+geom_histogram()+facet_grid(cl_annot_sex~.)
df_switch2 %>% filter(switching_category=="some_switch") %>%
  ggplot(aes(y=avg_p_male, x=frac_m))+geom_point(alpha=0.5)+
  #geom_errorbar(aes(ymin=ci_l, ymax=ci_u), alpha=0.5)+
  facet_grid(cl_annot_sex~.)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.5, lty=2, col="gray")+
  geom_hline(yintercept=1, lty=1, col="gray")+
  geom_hline(yintercept=0, lty=1, col="gray")



require('diptest')
install.packages('diptest')

dps <-samples_cl_sl %>% filter(!is.na(cl_acc)) %>% group_by(cl_acc) %>%
  mutate(dip_p=dip.test(p_male)$p)
summary(dps$dip_p)
dip.test(samples_cl_sl %>% filter(cl_name=="hff") %>% pull(p_male))

dps_some_swtch <- dps %>% filter(!is.na(cl_acc)) %>% select(cl_name, p_male, dip_p) %>% 
  right_join(df_switch2 
             %>% filter(frac_m > 0.40 & frac_m < 0.60) 
             ,
             by="cl_acc") %>%
  filter(switching_category=="some_switch") %>% 
  mutate(dip_pass=ifelse((dip_p > 0.05), "unimodal", "mulitmodal"))
table(dps_some_swtch %>% select(cl_acc, dip_pass) %>% unique() %>% pull(dip_pass))
76/(151+76)
# dip_pass == unimodal (null hypothesis is unimodal)

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

df_switch2 %>% left_join(cell_sex2 %>% select(cl_acc, cl_name)) %>%
  select(cl_acc, cl_name, cl_annot_sex, cl_allele_sex, avg_p_male,
         se_p_male, switching_category,
         num_s, frac_f, frac_m) %>%
  write_csv("data/cell_line_switching.csv")

# 98/2047
hc_switch <- df %>% 
  filter((cl_annot_sex=="female" & ci_l > 0.6) | 
           (cl_annot_sex=="male" & ci_u < 0.4),
         num_s > 3) %>%
  select(-expr_sex, -frac_f, -frac_m, -se_p_male) %>%
  mutate(novel=(cl_allele_sex == "unknown" |
                   cl_allele_sex==cl_annot_sex)) %>%
  arrange(avg_p_male) %>% 
  left_join(cell_sex2 %>% select(cl_acc, cl_name)) %>%
  mutate(cl_name=factor(cl_name, levels=cl_name)) 

table(hc_switch$novel, hc_switch$cl_annot_sex)
#female male
#FALSE      1   57
#TRUE       8   32

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

# try geom_label_repel?

### STOP ###

df %>% filter(annot_all %in% c("male>male", "male>unknown") &
                avg_p_male<0.5 & num_s > 1) %>%
  arrange(avg_p_male) %>% 
  mutate(cl_acc=factor(cl_acc, levels=cl_acc)) %>% # 30 
ggplot(aes(y=avg_p_male, x=cl_acc, col=annot_all))+
  geom_point(alpha=0.2, aes(size=num_s))+
  geom_errorbar(aes(ymin=(avg_p_male-1.96*se_p_male),
                    ymax=(avg_p_male+1.96*se_p_male)), 
                alpha=0.2,
                width=0.3)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("Average P(male)")+
  geom_hline(yintercept=0.5, lty=2, col="gray")+
  geom_hline(yintercept=1, lty=1, col="gray")+
  geom_hline(yintercept=0, lty=1, col="gray")



df %>% filter(annot_all=="male>unknown" &
                avg_p_male<0.5) # 98

df %>% filter(annot_all=="female>female" &
                avg_p_male>0.5) # 12

df %>% filter(annot_all=="female>unknown" &
                avg_p_male>0.5) # 31


df %>% filter(annot_all %in% c("female>female", "female>unknown") &
                avg_p_male>0.5 & num_s > 1) %>%
  arrange(desc(avg_p_male)) %>% 
  mutate(cl_acc=factor(cl_acc, levels=cl_acc)) %>%# 30 
  ggplot(aes(y=avg_p_male, x=cl_acc, col=annot_all))+
  geom_point(alpha=0.2, aes(size=num_s))+
  geom_errorbar(aes(ymin=(avg_p_male-1.96*se_p_male),
                    ymax=(avg_p_male+1.96*se_p_male)), 
                alpha=0.2,
                width=0.3)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("Average P(male)")+
  geom_hline(yintercept=0.5, lty=2, col="gray")+
  geom_hline(yintercept=1, lty=1, col="gray")+
  geom_hline(yintercept=0, lty=1, col="gray")

# ACK. I dont know how to best show this
f_hc_swap <- df %>% filter(annot_all %in% c("female>female", "female>unknown") &
                ((avg_p_male - 1.96*se_p_male) > 0.5 ) & num_s > 1)
f_hc_swap %>% arrange(desc(avg_p_male)) %>% View()
m_hc_swap <- df %>% filter(annot_all %in% c("male>male", "male>unknown") &
                             ((avg_p_male + 1.96*se_p_male) < 0.5 ) & num_s > 1)
m_hc_swap %>% arrange(avg_p_male) %>% View()

ggplot(cl_level_lab %>% 
         filter(!str_detect(cl_acc, ";")) %>%
         group_by(cl_acc) %>%
          mutate(annot_all=paste(cl_annot_sex,cl_allele_sex)), 
                               aes(x=cl_annot_sex, y=avg_p_male, 
                                   group=annot_all,
                                   col=cl_allele_sex))+
  geom_boxplot()
# problem: doesn't capture n per category
# ?? CONFUSION MATRIX??


### STOP ###


# // TODO: we want this on a *BY* cell line basis too!
# // TODO: what does multi-map mean if both are the same?
cl_data <- cl_df %>% select(cl_acc, cl_name, multi_map, cl_annot_sex, cl_allele_sex) %>% unique()
cl_data %>% filter(multi_map=="y") %>% unique()



head(multi_cl_acc_w_lab)

# counts?
counts_df <- cl_df %>% 
#  select(-organism, -p_male, -cl_acc, -metadata_sex, -study_acc, -cl_name) %>%
  mutate(cl_allele_sex=ifelse(multi_map=="y", "multi_map", cl_allele_sex),
         cl_annot_sex=ifelse(multi_map=="y", "multi_map", cl_annot_sex)) %>%
  group_by(data_type, cl_annot_sex,cl_allele_sex, expr_sex) %>%
  count() %>%
  ungroup() %>%
  replace_na(list("expr_sex"="unknown"))


present_data <- counts_df %>% 
  filter(cl_annot_sex != "multi_map" & 
                       cl_annot_sex != "unknown" &
                       expr_sex != "unknown")
(num_switching <- present_data %>% 
  filter(cl_annot_sex=="male" &
           expr_sex=="female") %>%
  summarize(sum(n)))

switching_prev <- present_data %>% 
  filter(cl_annot_sex=="male" &
           expr_sex=="female" & (expr_sex==cl_allele_sex | cl_allele_sex=="both"))
switching_new <- present_data %>%
  filter(cl_annot_sex=="male" &
           expr_sex=="female" & (expr_sex!=cl_allele_sex & cl_allele_sex!="both"))
(num_switching_prev <- switching_prev %>%
  summarize(sum(n)))

num_switching_prev/num_switching
# this is number of samples :/ not number of cell lines


### TOOO SLOWWWWW
# require(tictoc) # should take like 7 mins. agh
# tic()
# cl_multi <- cl %>% 
#   semi_join(cl2 %>% filter(n>1)%>% sample_n(400)) %>%
#   group_by(sample_acc, organism, data_type, source_type, metadata_sex, expr_sex) %>%
#   select(-p_male, -study_acc) %>%
# #  summarize_at(vars(cl_acc, cl_name, cl_annot_sex, cl_allele_sex), 
# #               paste(unique(.), collapse=";"))
#   mutate_at(vars(cl_acc, cl_name, cl_annot_sex, cl_allele_sex), 
#             ~paste(unique(.), collapse=";"))
#   #       cl_name=paste(unique(cl_name), collapse=";"),
#   #      cl_annot_sex=paste(unique(cl_annot_sex), collapse=";"),
#   #      cl_allele_sex=paste(unique(cl_allele_sex), collapse=";")) 
# toc()
# 
# 
# stopifnot(nrow(samples_cell_sl2)==nrow(mu_s3))
# samples_cell_sl2 %>% write_csv("data/cl_sample_expr_sex.csv")      
# 
# # --- 
