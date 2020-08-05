# Code for getting a list of sample source types using the metadata

# ----- grab cell labeling info ----- #
refine_mapped <- read_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl.csv", 
                                  prefix, data_type2))
refine_mapped2 <- read_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl_part2.csv",
                                   prefix, data_type2))

cl_stud <- fread(sprintf("data/02_labeled_data/%s_%s_study_cl.csv", prefix, data_type2), data.table = FALSE)

no_cell <- metadata %>% filter(cl_line=="" | cl_line=="--")  %>% select(acc)
some_cl <- metadata %>% filter(cl_line!="" | cl_line!="--")  

cl_lab_unique <- refine_mapped  %>% 
  bind_rows(refine_mapped2 %>% select(gsm, orig_str, accession)) %>% 
  unique()

xenograft_cl <- cl_lab_unique %>% filter(str_detect(orig_str, "xenograft"))
primary_cl <- cl_lab_unique %>% filter(str_detect(orig_str, "primary"))
stem_cl <- cl_lab_unique %>% filter(str_detect(orig_str, "stem|ipsc"))

meta_unique <- metadata %>% select(acc, cl_line, part, title) 

#unmapped_cl <- some_cl %>% anti_join(cl_lab_unique, by=c("acc"="gsm"))

mu_s <- meta_unique %>% left_join(cl_lab_unique, by=c("acc"="gsm")) %>% 
  group_by(acc)  %>% mutate(str=paste(c(cl_line, part, title), collapse=" ")) %>%
  ungroup() %>%
  mutate(str=tolower(str)) # SLOW
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

# now add in the study_level
samp_w_study_cl <- study_cl %>% select(gse, accession) %>% unique() %>% left_join(mapping, by=c("gse"="study_acc")) 
table((mu_s2 %>% semi_join(samp_w_study_cl, by=c("acc"="sample_acc")))$source_type)

# ok most of these already assigned
mu_s3 <- mu_s2 %>% mutate(source_type=ifelse(source_type=="tissue" & acc %in% samp_w_study_cl$sample_acc, "other", source_type))
table((mu_s3 %>% semi_join(samp_w_study_cl, by=c("acc"="sample_acc")))$source_type)
sample_source_df <- mu_s3 %>% select(acc, source_type)
sample_source_df %>% write_csv(sprintf("data/02_labeled_data/%s_%s_sample_source.csv", prefix, data_type2))
