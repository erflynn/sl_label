require('data.table')
cell_lab <- fread("../labeling/geo2drug/data/00_db_data/cellosaurus_df.txt", data.table=FALSE)
cl_sex <- cell_lab %>% mutate(accession=tolower(accession), annot_sex=tolower(sex)) %>% select(accession, annot_sex) 
h_cl_sex <- h_cell %>% filter(!is.na(accession)) %>% separate_rows(accession, sep=",") %>% left_join(cl_sex) %>%
  mutate(annot_sex=ifelse(annot_sex == "" | annot_sex=="sex unspecified", "no_annot", annot_sex))
table(h_cl_sex$annot_sex, h_cl_sex$expr_sex)
sum(h_cl_sex$annot_sex== h_cl_sex$expr_sex)/nrow(h_cl_sex %>% filter(annot_sex!="no_annot"))
counts_cl_sex <- h_cl_sex %>% 
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