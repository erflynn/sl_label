require('data.table')
cell_lab <- fread("data/00_db_data/cellosaurus_df_v2.txt", data.table=FALSE)
cl_sex <- cell_lab %>% 
  mutate(accession=tolower(primary_accession), annot_sex=tolower(sex), alleles) %>% 
  select(accession, annot_sex, alleles) 
refine_mapped <- read_csv(sprintf("data/%s_acc_to_cl.csv", prefix))

human_sl <- read_csv("data/02_labeled_data/human_all_sl.csv")
human_metadata <- read.csv("data/01_metadata/human_metadata.csv", stringsAsFactors = FALSE)

h_cl_sex <- refine_mapped %>% filter(!is.na(accession)) %>% separate_rows(accession, sep=",") %>% 
  inner_join(cl_sex) %>%
  mutate(annot_sex=ifelse(annot_sex == "" | annot_sex=="sex unspecified", 
                          "no_annot", annot_sex))
h_cl_sex2 <- h_cl_sex %>% left_join(human_sl, by=c("gsm"="id")) %>% 
  mutate(expr_sex=ifelse(pred>0.5, "male", "female"))
table(h_cl_sex2$annot_sex, h_cl_sex2$expr_sex)
sum(h_cl_sex2$annot_sex== h_cl_sex2$expr_sex)/nrow(h_cl_sex2 %>% filter(annot_sex!="no_annot"))
mismatch <- h_cl_sex2 %>% filter(annot_sex != expr_sex & annot_sex !="no_annot") %>%
  mutate(allele_sex=case_when(
    alleles=="" ~ "",
    alleles=="X" ~ "female",
    alleles=="X,Y" ~ "male",
    alleles %in% c("X|X,Y", "X,Y|X") ~ "both",
    alleles == "Not_detected|X" ~ "female",
    TRUE ~ alleles
  ))
mismatch_new <- mismatch %>% filter(allele_sex!=expr_sex & allele_sex!="both")
hc_mismatch_new <- mismatch_new %>% filter(pred < 0.2 | pred > 0.8) # 155 cell lines, 3717 mismatches
hc_mismatch_new 

# add in the study
# is the study mixed sex?
exp_sample <- read_csv("data/01_metadata/human_exp_to_sample.csv")
hc_mismatch_new2 <- hc_mismatch_new %>% left_join(exp_sample, by=c("gsm"="sample_acc")) %>% rename(gse=study_acc) %>%
  select(gse, gsm, everything())

load("data/summary_files.RData")

hc_2 <- hc_mismatch_new2 %>% left_join(h_summ %>% filter(labeling_method=="exprsex") %>% select(study, sex) %>%
                                        rename(gse=study, study_sex=sex) %>% mutate(gse=as.character(gse)), by=c("gse"))
hc_2 %>% write_csv("data/mismatch_new_cl.csv")
hc_2 %>% filter(study_sex=="mixed" & expr_sex=="male")

hc_3 <- hc_2 %>% select(gse) %>% unique() %>% left_join(exp_sample, by=c("gse"="study_acc")) %>%
  left_join(human_metadata %>% select(acc, cl_line), by=c("sample_acc"="acc"))  %>%
  rename(gsm=sample_acc)
hc_3 %>% group_by(gse) %>% 
  summarize(num_samples=n(), cl_line=paste(setdiff(unique(cl_line[!is.na(cl_line)]), c("", "NA")), collapse="|")) %>%
  left_join(h_summ %>% filter(labeling_method=="exprsex") %>% select(study, sex) %>% rename(study_sex=sex), by=c("gse"="study")) %>%
  write_csv("data/study_cl_line.csv")

hc_3 %>% select(gse, gsm) %>% left_join(h_cl_sex2, by=c("gsm")) %>%
  write_csv("data/study_sample_level_cl.csv")

counts_cl_sex <- h_cl_sex2 %>% 
  filter(annot_sex != "no_annot") %>% 
  select(gse, accession, annot_sex, expr_sex) %>%
  group_by(gse, accession, annot_sex) %>%
  summarize(num_samples=n(),  
            num_correct=sum(expr_sex==annot_sex), num_incorrect=sum(expr_sex!=annot_sex)) %>%
  mutate(frac_incorrect=num_incorrect/num_samples) %>%
  arrange(desc(frac_incorrect))

counts_cl_sex %>% arrange(desc(num_incorrect))

table(counts_cl_sex$frac_incorrect==1) # 107
table(counts_cl_sex$frac_incorrect==0) # 327

counts_cl_sex %>% filter(frac_incorrect==1) %>% arrange(desc(num_incorrect))


cell_annot %>% 
  filter(gse=="GSE27211")

cell_annot %>% 
  filter(gse=="GSE73935")

h_cell %>% filter(gse=="GSE27211") %>% select(gsm, expr_sex) %>% left_join(human, by=c("gsm"="acc")) %>% head()

# // TODO - cell_line isn't correct! need a combo :(



refine_cl_annot <- human_metadata %>% select(acc, cl_line) %>% filter(!is.na(cl_line) & cl_line!="--" & cl_line !="")

counts_cl_sex %>% filter(accession=="cvcl_0302") %>% arrange(desc(num_samples)) %>% write_csv("data/cvcl_0302.csv") # 23 studies

counts_cl_sex %>% filter(accession=="cvcl_0302") %>% ungroup() %>% group_by(accession, annot_sex) %>% summarize(num_studies=n(), num_samples=sum(num_samples), num_incorrect=sum(num_incorrect))


### next step -- with cell line annot
res2 <- read_csv("cell_line_no_trt.csv")
human_metadata <- read.csv("data/human_metadata2.csv")
res2.2 <- res2 %>% left_join(human_metadata %>% select(acc, sex_lab)) 
res2.3 <- res2.2 %>% rename(expr_sex=sex_lab) %>% left_join(cl_sex) %>% mutate(incorr=(annot_sex!=expr_sex))
table(res2.3$annot_sex==res2.3$expr_sex) # 8981 true, 1020 false

res2.3 %>%  group_by(accession, incorr) %>% count() %>% arrange(desc(n)) %>% ungroup() %>% pivot_wider(id_cols=accession, names_from=incorr, values_from=n) %>% rename("num_corr"="FALSE", "num_inc"="TRUE") %>% mutate(frac_inc=num_inc/(num_corr+num_inc)) %>% arrange(desc(frac_inc))
# HT116 and A549, LNCAP