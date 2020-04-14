# divide into cell line, primary cells, stem cells, iPSCs, and tissues

require('tidyverse')

sample_metadata <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", prefix))
rnaseq_sample_metadata <- read.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", prefix))

# put these together
meta_unique <- sample_metadata %>% select(acc, cl_line, part, title) %>% bind_rows(
  rnaseq_sample_metadata %>% select(acc, cl_line, part, title)) %>%
  unique() # 560297

compendia_cl <- read_csv("data/02_labeled_data/human_compendia_sample_cl.csv")
rnaseq_cl <- read_csv("data/02_labeled_data/human_rnaseq_sample_cl.csv")
cl_lab_unique <- bind_rows(compendia_cl, rnaseq_cl) %>% unique() # 63843
cl_lab_unique %>% head()
xenograft_cl <- cl_lab_unique %>% filter(str_detect(orig_str, "xenograft"))
primary_cl <- cl_lab_unique %>% filter(str_detect(orig_str, "primary"))
stem_cl <- cl_lab_unique %>% filter(str_detect(orig_str, "stem|ipsc"))

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

# now add the sex labels and look at the switching
compendia_sex <- read.csv("data/01_metadata/human_microarray_metadata_sex.csv") %>% unique()
rnaseq_sex <- read.csv("data/01_metadata/human_rnaseq_metadata_sex.csv") %>% unique()
metadata_sex <- compendia_sex %>% bind_rows(rnaseq_sex) %>% unique()

mu_s3 <- mu_s2 %>% 
  left_join(metadata_sex %>% select(acc, mapped_sex), by=c("acc")) %>% 
  rename(metadata_sex=mapped_sex)

# for now we use the compendia labels unless they are missing
compendia_sl <- read_csv("data/02_labeled_data/human_all_sl.csv") 
rnaseq_sl <- read_csv("data/02_labeled_data/human_rnaseq_sl.csv") 
all_sl <- compendia_sl %>% bind_rows(rnaseq_sl %>% anti_join(compendia_sl, by=c("id")) ) %>% 
  mutate(exprsex=ifelse(pred > 0.5, "male", "female"))

mu_s4 <- mu_s3 %>% left_join(all_sl, by=c("acc"="id"))
tissue <- mu_s4 %>% filter(source_type=="tissue")
primary <- mu_s4 %>% filter(source_type=="primary_cells")
named_cl <- mu_s4 %>% filter(source_type=="named_cl")
# also grab the annot_sex for the named_cl
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

# add the annotated sex
named_cl_df <- named_cl %>% select(acc, accession, metadata_sex, pred, exprsex) %>% 
  separate_rows(accession, sep=";") %>%
  left_join(cell_sex, by=c("accession"="primary_accession")) 

named_cl_df2 <- named_cl_df %>%
  group_by(acc) %>%
  summarise_all(~paste(sort(unique(unlist(.))), collapse=";"))
named_cl_df3 <- named_cl_df2 %>% filter(!is.na(pred)) 

# is this ever different from the metadata sex??
mismatch <- named_cl_df3 %>% filter(annot_sex != metadata_sex & metadata_sex != "unknown") # 238
mismatch %>% View()

# ones that match
matched <- named_cl_df3 %>% filter(annot_sex == metadata_sex | metadata_sex== "unknown")
matched2 <- matched %>% filter(exprsex!="") %>% 
  mutate(annot_sex2=ifelse(annot_sex %in% c("male", "female"), annot_sex, "unknown")) %>%
  mutate(source_type="cell_line") %>%
  select(-metadata_sex) %>%
  rename(metadata_sex=annot_sex2) 
table(matched2$metadata_sex, matched2$exprsex)

primary2 <- primary %>% select(acc, metadata_sex, pred, exprsex, source_type) 



table(primary2$metadata_sex, primary2$exprsex)
chisq.test(matrix(c(452, 17, 62, 623), nrow=2,byrow = TRUE))
tissue2 <- tissue %>% select(acc, metadata_sex, pred, exprsex, source_type) 
table(tissue2$metadata_sex, tissue2$exprsex)

df <- do.call(rbind, list(
matched2 %>% select(acc, metadata_sex, exprsex, pred, source_type),
primary2 %>% select(acc, metadata_sex, exprsex, pred, source_type),
tissue2 %>% select(acc, metadata_sex, exprsex, pred, source_type))) %>% filter(!is.na(source_type)) %>%
  filter(!is.na(exprsex))
df2 <- df %>% pivot_longer(c(metadata_sex, exprsex), names_to="labeling_method", values_to="sex") %>%
  mutate(labeling_method=ifelse(labeling_method=="metadata_sex", "metadata", "expression"))

require("ggalluvial")
alluvPlotSample <- function(ds){
  flow_freq_counts <- ds %>% 
    ungroup() %>%
    rename(metadata=metadata_sex, expression=exprsex) %>%
    group_by(metadata, expression) %>% 
    mutate(Freq=n()) %>% 
    select(-acc, -source_type, -pred) %>% 
    unique() %>%
    ungroup() %>%
    mutate(row_id=1:n()) %>%
    gather(key="labeling_method", value="sex", -Freq, -row_id) %>%
    mutate(row_id=as.factor(row_id), 
           labeling_method=factor(labeling_method, levels=c("metadata", "expression")),
           sex=as.factor(sex)) %>%
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
    theme(legend.position = "none") 
}
alluvPlotSample(df %>% filter(source_type=="cell_line" ))
alluvPlotSample(df %>% filter(source_type=="primary_cells" ))
alluvPlotSample(df %>% filter(source_type=="tissue" & metadata_sex!="mixed" ))

df_label <- df %>% filter(!metadata_sex %in% c("unknown", "mixed") & metadata_sex!=exprsex) %>% 
  mutate(label=paste(metadata_sex, exprsex, sep=">"))
ggplot(df_label, aes(x=source_type, fill=label))+geom_bar(position="fill")+ylab("fraction")

matched2 %>% head()

# distribution of pred scores?
df3 <- df %>% mutate(pred=as.numeric(as.character(pred))) #%>% filter(metadata_sex!=exprsex) 

ggplot(df3, aes(x=pred))+geom_density(aes(col=source_type))+xlab("P(male)")

# look at how the mislabeled ones pan out


# alluvial with allele sex"
matched3 <- matched2 %>% mutate(allele_sex=case_when(
  (allele_sex %in% c("female")) ~ "female",
  (allele_sex %in% c("male")) ~ "male",
  allele_sex=="both" ~ "both",
  allele_sex=="unknown" ~ "unknown",
  TRUE ~ "multi-map"
  )) %>%
  filter(!is.na(exprsex) & !is.na(allele_sex) & !is.na(metadata_sex)) %>%
  select(-accession, -alleles, -annot_sex, -cl)


flow_freq_counts <- matched3 %>% filter(source_type=="cell_line") %>%
  ungroup() %>%
  rename(metadata=metadata_sex, expression=exprsex, amelogenin=allele_sex) %>%
  group_by(metadata, expression, amelogenin) %>% 
  mutate(Freq=n()) %>% 
  select(-acc, -source_type, -pred) %>% 
  unique() %>%
  ungroup() %>%
  mutate(row_id=1:n()) %>%
  gather(key="labeling_method", value="sex", -Freq, -row_id) %>%
  mutate(row_id=as.factor(row_id), 
         labeling_method=factor(labeling_method, levels=c("metadata", "amelogenin", "expression")),
         sex=as.factor(sex)) %>%
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
  theme(legend.position = "none") 

# what fraction of m--> f have amelogenin
matched4 <- matched3 %>% filter(!allele_sex %in% c("unknown", "multi-map") & !metadata_sex=="unknown")
matched4 %>% filter(metadata_sex=="male" & exprsex=="female") %>% nrow() # 8417
matched4 %>% filter(metadata_sex=="male" & exprsex=="female" & 
                      allele_sex %in%c("female","both")) %>% nrow() # 6625

matched4 %>% filter(metadata_sex=="male" & allele_sex %in%c("female","both")) %>% nrow() # 8703
matched4 %>% filter(metadata_sex=="female" & allele_sex %in%c("male","both")) %>% nrow() # 53

matched4 %>% filter(metadata_sex=="female" & exprsex=="male") %>% nrow() # 1015
matched4 %>% filter(metadata_sex=="female" & exprsex=="male" & 
                      allele_sex %in%c("male","both")) %>% nrow() # 12


### look at annie's ones
exp_to_samp <- read_csv("data/01_metadata/human_exp_to_sample.csv")
exp_to_samp_rnaseq <- read_csv("data/01_metadata/human_exp_to_sample_counts.csv")
all_to_study <- exp_to_samp %>% bind_rows(exp_to_samp_rnaseq %>% select(study_acc, sample_acc)) %>% unique()
#study1 <- exp_to_samp %>% filter(study_acc=="SRP095272") %>% 
# left_join(matched4, by=c("sample_acc"="acc")) %>%
#  filter(!is.na(pred))

studies <- exp_to_samp %>% filter(study_acc%in% c(
  "SRP095272","GSE48398", "GSE25429", "GSE47856", "GSE109021", "SRP134389") )%>% 
  left_join(named_cl_df3, by=c("sample_acc"="acc")) %>%
  filter(!is.na(pred)) %>%
  mutate(pred=as.numeric(pred)) 
ggplot(studies %>% filter(annot_sex %in% c("female", "male")), 
       aes(x=pred, y=annot_sex))+
  geom_point(aes(col=annot_sex), alpha=0.3, position=position_jitter())+
  facet_wrap(.~study_acc, ncol=1)+
  xlab("P(male)")

studies %>% filter(study_acc=="SRP095272") %>% filter(exprsex!=annot_sex) %>% View()

studies %>% filter(study_acc=="SRP134389") %>% filter(exprsex!=annot_sex) %>% View()
studies %>% filter(study_acc=="SRP134389") %>% select(cl) %>% unique()


#####
t2 <- tissue2 %>% filter(metadata_sex %in% c("male", "female") & !is.na(pred)) %>% 
  left_join(all_to_study, by=c("acc"="sample_acc")) 


