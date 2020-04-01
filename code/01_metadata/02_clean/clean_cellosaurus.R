require('tidyverse')
require('data.table')

# look at this: http://metasra.biostat.wisc.edu/static/publication_datasets/cvcl_mappings.json
# ATCC cell lines?


cell_info_df <- fread("data/00_db_data/cellosaurus_df_v2.txt", data.table=FALSE)

cell_lab <- cell_info_df %>% 
  select(synonyms, primary_accession, cl) %>%
  rename(accession=primary_accession)
cell_lab2 <- data.frame(apply(cell_lab, c(1,2) , function(x) str_trim(tolower(x))))
cell_lab3 <- cell_lab2 %>% separate_rows(synonyms, sep="\\|")  %>%
  mutate(synonyms=ifelse(synonyms==cl, "", synonyms))
cell_lab_syn <- cell_lab3 %>% select(accession, synonyms) %>% 
  filter(synonyms != "")  %>% unique()
cell_lab_name <- cell_lab3 %>% select(accession, cl) %>% distinct()

duplicated_names <- cell_lab_name$cl[duplicated(cell_lab_name$cl)] # 31
duplicated_synonyms <- cell_lab_syn$synonyms[duplicated(cell_lab_syn$synonyms)] # 760 
dup_name_syn <- intersect(cell_lab_name$cl, cell_lab_syn$synonyms) # 412

cell_lab_name2 <- cell_lab_name %>% 
  unique() %>%
  mutate(new_cl=str_replace_all(cl, "[-|\\.|/| |#|_|(|)|+]", ""))

duplicated_punct_diff <- cell_lab_name2$new_cl[duplicated(cell_lab_name2$new_cl)] # 236
dup_punct_diff2 <- setdiff(duplicated_punct_diff, duplicated_names) # 294

cell_lab_syn_alt <- cell_lab_syn %>% 
  mutate(new_syn=str_replace_all(synonyms, "[-|\\.|/| |#|_|(|)|+]", "")) %>%
  select(accession, new_syn) %>%
  unique()

dup_punct_diff_syn <- cell_lab_syn_alt$new_syn[duplicated(cell_lab_syn_alt$new_syn)]  # 913
dup_punct_syn2 <- setdiff(dup_punct_diff_syn, duplicated_synonyms) # 579

# check for overlap
name_syn <- cell_lab_name2 %>%
  left_join(cell_lab_syn_alt, by=c("accession")) %>%
  filter(new_cl != new_syn & !is.na(new_syn))
multi_syn_name <- intersect(unique(name_syn$new_cl), unique(name_syn$new_syn))
dup_punct_syn_name <- setdiff(multi_syn_name, dup_name_syn) # 91

# synonyms and multi-mapping
name_df <- cell_lab_name2 %>% select(cl, accession) %>% unique() %>% 
  bind_rows(cell_lab_name2 %>% select(new_cl, accession) %>% unique() %>% rename(cl=new_cl)) %>%
  unique() %>%
  group_by(cl) %>%
  summarize(accession=paste(accession, collapse="|")) # 159735

# remove all that map to a name
syn_df <- cell_lab_syn_alt %>% rename(cl=new_syn) %>%
  bind_rows(cell_lab_syn %>% rename(cl=synonyms)) %>%
  unique() %>%
  group_by(cl) %>%
  summarize(accession=paste(accession, collapse="|")) # 92606

syn_df_filt <- syn_df %>% anti_join(name_df, by=c("cl")) # 85075

mult_name <- name_df %>% filter(str_detect(accession, "\\|"))
mult_syn <- syn_df_filt %>% filter(str_detect(accession, "\\|"))

mult_name2 <- mult_name %>% rename(acc=accession) %>%
  separate_rows(acc, sep="\\|") %>%
  inner_join(cell_info_df %>% 
               select(primary_accession, derivs, origins, accession) %>%
               mutate_all(tolower), 
             by=c("acc"="primary_accession")) %>%
  separate_rows(accession, sep="\\|") %>%
  filter(accession != "list()") %>%
  separate_rows(derivs, origins, sep="\\|") %>%
  group_by(cl, acc) %>%
  summarize(addl_acc=paste(setdiff(unique(unlist(c(derivs, origins, accession))), c("")), collapse="|"))


mult_syn2 <- mult_syn %>% rename(acc=accession) %>%
  separate_rows(acc, sep="\\|") %>%
  inner_join(cell_info_df %>% 
               select(primary_accession, derivs, origins, accession) %>%
               mutate_all(tolower), 
             by=c("acc"="primary_accession")) %>%
  separate_rows(accession, sep="\\|") %>%
  filter(accession != "list()") %>%
  separate_rows(derivs, origins, sep="\\|") %>%
  group_by(cl, acc) %>%
  summarize(addl_acc=paste(setdiff(unique(unlist(c(derivs, origins, accession))), c("")), collapse="|"))


mult_name3 <- mult_name2 %>% 
  group_by(cl) %>% 
  mutate(overlap=paste(unlist(strsplit(addl_acc, "\\|"))[duplicated(unlist(strsplit(addl_acc, "\\|")))], 
                       collapse="|" ))  %>% 
  mutate(res_ids=ifelse(overlap=="", acc, overlap))

mult_name4 <- mult_name3 %>% 
  group_by(cl) %>%
  summarize(list_acc=paste(acc, collapse="|"), 
            list_dedup_acc=paste(unique(res_ids), collapse="|"))

mult_syn3 <- mult_syn2 %>%
  group_by(cl) %>% 
  mutate(overlap=paste(unlist(strsplit(addl_acc, "\\|"))[duplicated(unlist(strsplit(addl_acc, "\\|")))], 
                       collapse="|" ))  %>% 
  mutate(res_ids=ifelse(overlap=="", acc, overlap))

mult_syn4 <- mult_syn3 %>% 
  group_by(cl) %>%
  summarize(list_acc=paste(acc, collapse="|"), 
            list_dedup_acc=paste(unique(res_ids), collapse="|"))

# how many map to parents
syn_parent <- mult_syn3 %>% filter(overlap!="" & acc != addl_acc) 
name_parent <- mult_name3 %>% filter(overlap!="" & acc != addl_acc) 
name_parent_df <- mult_name4 %>% semi_join(name_parent) # 77
syn_parent_df <- mult_syn4 %>% semi_join(syn_parent) # 739

mult_name5 <- mult_name4 %>% group_by(cl) %>% mutate(
  nwords=length(strsplit(cl, " ")[[1]]),
  numchar=nchar(cl))

mult_syn5 <- mult_syn4 %>% group_by(cl) %>% mutate(
  nwords=length(strsplit(cl, " ")[[1]]),
  numchar=nchar(cl))
  
mult_name5 %>% write_csv("data/00_db_data/cell_line/name_to_acc.csv")
mult_syn5 %>% write_csv("data/00_db_data/cell_line/syn_to_acc.csv")

name_parent_df %>% write_csv("data/00_db_data/cell_line/name_parent.csv")
syn_parent_df %>% write_csv("data/00_db_data/cell_line/syn_parent.csv")

mult_cell_df <- rbind(mult_name5 %>% mutate(type="name"), mult_syn5 %>% mutate(type="syn")) %>% 
  select(cl, list_dedup_acc, nwords, numchar, type) %>%
  rename(acc=list_dedup_acc) # 1321

# now add back in the original stuff
cell_df1 <- name_df %>% mutate(type="name") %>%
  bind_rows(syn_df_filt %>% mutate(type="syn")) %>%
  anti_join(mult_cell_df, by=c("cl")) %>%
  rename(acc=accession) %>%
  mutate(
    nwords=length(strsplit(cl, " ")[[1]]),
    numchar=nchar(cl)) 
cell_df <- bind_rows(mult_cell_df, cell_df1)  

cell_df %>% write_csv("data/00_db_data/cell_syn_df.csv")
