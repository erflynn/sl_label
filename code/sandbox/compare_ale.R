ale_sex_lab <- read_csv(sprintf("../sex_labeling/geo_pipeline/data/01_sample_lists/%s_ale_sex_lab.csv", prefix))
feat.df3 <- feat.df3 %>% mutate(acc=as.character(acc), sex=as.character(sex))

length(intersect(ale_sex_lab$gsm, feat.df3$acc)) # 50,042 and 116,318 human
# mouse: there are 176,035 ale, and 279,998 in refinebio... hmm
# maybe this is because of differences in sample ids?
# -ALSO- all the ALE ones *HAVE* labels
ale2 <- ale_sex_lab %>% filter(gsm %in%  feat.df3$acc )
ale3 <- ale2 %>% left_join(feat.df3 %>% select(acc, sex), by=c("gsm"="acc"))
ale4 <- ale3 %>% mutate(ale_sex=case_when(
  Gender=="M" ~ "male",
  Gender=="F" ~ "female",
  TRUE ~ Gender)) %>%
  rename(rb_sex=sex)
ale4 %>% filter(rb_sex==ale_sex & rb_sex!="") %>% nrow()
# mouse: 37151, 74.2% agreement...
# human: 92094, 79.2% agreement

ale4 %>% filter(rb_sex!=ale_sex) %>%
  select(rb_sex, ale_sex) %>% unique()
# mouse: 78 diff mismatches
# human: 438 diff mismatches

# BUT - there are only two exact swaps!
ale4 %>% filter(ale_sex=="male" & rb_sex=="female")
ale4 %>% filter(ale_sex=="female" & rb_sex=="male")
# GSE99192 - this is the study, and rb is wrong... (superseries + subseries)
# GSE17743, GSE55232, GSE83951 -- all of them rb says male and ale says female

# -----  create a training list ----- #

# label columns by file + filter for cols actually present
m_cols <- fread(sprintf("data/%s/02_sample_lists/%s_cols.txt",prefix, prefix),
                data.table=FALSE, stringsAsFactors = FALSE)
m_to_idx <- data.frame(t(m_cols))
colnames(m_to_idx) <- c("col_name")
rownames(m_to_idx) <- NULL
m_to_idx <- m_to_idx %>%
  mutate(idx=1:nrow(m_to_idx)) %>%
  mutate(f_idx=(idx-2)%/%CHUNK.SIZE) %>%
  filter(!is.na(col_name))

sex_lab2 <- sex_lab %>% inner_join(m_to_idx, by=c("acc"="col_name"))
sex_lab2 %>%
  write_csv(sprintf("data/%s/02_sample_lists/%s_sex_lab_no_cl_idx.csv", prefix, prefix))

# randomly sample 1000 males and 1000 females
set.seed(228)
f_sample <- sex_lab2 %>%
  filter(sex=="female") %>%
  sample_n(1000)
m_sample <- sex_lab2 %>%
  filter(sex=="male") %>%
  sample_n(1000)
f_sample %>% write_csv(sprintf("data/%s/02_sample_lists/%s_f_sample.csv", prefix, prefix))
m_sample %>% write_csv(sprintf("data/%s/02_sample_lists/%s_m_sample.csv", prefix, prefix))