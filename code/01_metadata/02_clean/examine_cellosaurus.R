require('tidyverse')
require('data.table')
options(stringsAsFactors = FALSE)

# look at the breakdown of

# read in the list of repeated cell lines
t1 <- read.delim("data/00_db_data/cellosaurus_identical_names_t1.txt", sep=":", header=FALSE)
t2 <- read.delim("data/00_db_data/cellosaurus_identical_names_t2.txt", sep=":", header=FALSE)
t3 <- read.delim("data/00_db_data/cellosaurus_identical_names_t3.txt", sep=":", header=FALSE)
t4 <- read.delim("data/00_db_data/cellosaurus_identical_names_t4.txt", sep=":", header=FALSE)

t1.2 <- t1 %>% rename(cl=V1, acc=V2) %>% 
  separate_rows(acc, sep=",") %>% 
  mutate(acc=str_trim(acc)) %>%
  mutate(cl=tolower(cl))
t2.2 <- t2 %>% rename(cl=V1, acc=V2) %>% 
  separate_rows(acc, sep=",") %>% 
  mutate(acc=str_trim(acc)) %>%
  group_by(acc) %>%
  mutate(cl=tolower(str_split(cl, pattern=",")[[1]][[1]])) %>%
  ungroup()
t3.2 <- t3 %>% rename(cl=V1, acc=V2) %>% 
  separate_rows(acc, sep=",") %>% 
  mutate(acc=str_trim(acc)) %>%
  group_by(acc) %>%
  mutate(cl=tolower(str_split(cl, pattern=",")[[1]][[1]])) %>%
  ungroup()
t4.2 <- t4 %>% rename(cl=V1, acc=V3, synonym=V4) %>% 
  select(-V2) %>%
  mutate(acc=str_replace_all(acc, ", Synonym", "")) %>%
  mutate(acc=str_trim(acc), synonym=str_trim(synonym)) %>%
  mutate(cl=tolower(cl))

# put 1-3 together

t_all <- do.call(rbind, list(t1.2, t2.2, t3.2)) %>% arrange(cl) %>% unique()
# 1063 (1061) accessions map to 500 cell names

cell_info_df <- fread("data/00_db_data/cellosaurus_df_v2.txt", data.table=FALSE)
derivs <- cell_info_df %>% filter(derivs != "" | origins != "")

# possibilities:
#  - same parent or origin, shared accession
#  - one is parent/origin of the other, or v.v.
#  - they're connected somehow on an insane tree
#  - no relationship


# now see if we can separate out the issues
t_all %>% head()
# t1_w_info <- t_all %>% inner_join(cell_info_df, by=c("acc"="primary_accession"))
# t1_w_info$acc2 <- sapply(t1_w_info$accession, function(x) strsplit(x, "\\|")[[1]][[2]])
# mult_acc <- t1_w_info %>% filter(acc2 != "list()") %>% select(cl.x, accession, acc, derivs, origins, acc2)
# # there are 9 with multiple accessions
# 
# t1_w_info2 <- t1_w_info %>% filter((derivs!="" & !is.na(derivs))) %>%
#   separate_rows(derivs, sep="\\|")
# 
# # now only keep the ones that are more than 1 deriv per cl.x, acc
# t1_w_info2.1 <- t1_w_info %>% semi_join(t1_w_info2 %>% 
#                                           select(cl.x, acc) %>% 
#                                           unique() %>% 
#                                           group_by(cl.x) %>% 
#                                           count() %>% 
#                                           filter(n > 1))
# 
# t1_w_info2.2 <- t1_w_info %>% filter((origins!="" & !is.na(origins))) %>%
#   separate_rows(origins, sep="\\|")
# t1_w_info2.21 <- t1_w_info2.2 %>% semi_join(t1_w_info2.2 %>% select(cl.x, acc) %>% 
#                                               unique() %>% 
#                                               group_by(cl.x) %>% 
#                                               count() %>% 
#                                               filter(n > 1))
# 
# t1_w_info3 <- t1_w_info2.1 %>% select(acc, accession, derivs, cl.x, cl.y) %>% arrange(cl.x) 
# t1_w_info3.2 <- t1_w_info2.21 %>% select(acc, accession, origins, cl.x, cl.y) %>% arrange(cl.x) %>% unique()
# t1_w_info4 <- t1_w_info3 %>% group_by(cl.x) %>% mutate(num_deriv=n()) %>% 
#   select(cl.x, derivs, num_deriv) %>% unique() %>% mutate(num_unique=n())
# t1_w_info4.2 <- t1_w_info3.2 %>% group_by(cl.x) %>% mutate(num_origins=n()) %>% 
#   select(cl.x, origins, num_origins) %>% unique() %>% mutate(num_unique=n())
# 
# same_parent <- t1_w_info4 %>% filter(num_deriv > num_unique) # these have the same parent!!
# diff_parent <-t1_w_info4 %>% filter(num_deriv == num_unique & num_unique != 1)
# t1_w_info4.2 %>% filter(num_origins > num_unique) # none of these
# 
# 
# same_parent_df <- t1_w_info3 %>% semi_join(same_parent, by=c("cl.x", "derivs"))
# length(unique(same_parent_df$derivs)) # 81 --> 48 cl names and 26 derivs (appears many are derivs)
# 
# # is any the parent of the other?
# 
# df_parent <- do.call(rbind,
#                      list(t1_w_info %>% select(cl.x, acc) %>% mutate(new_acc=acc), 
#                           t1_w_info %>% select(cl.x, acc, derivs) %>% filter(derivs != "") %>% rename(new_acc=derivs),
#                           t1_w_info %>% select(cl.x, acc, origins) %>% filter(origins != "")%>% rename(new_acc=origins))) %>% 
#   separate_rows(new_acc, sep="\\|")  %>%
#   unique()
# 
# df_parent %>% arrange(cl.x) %>% head()
# df_parent2 <- df_parent %>% unique() %>% group_by(cl.x) %>% mutate(nr=n()) %>% 
#   select(cl.x, new_acc, nr) %>% unique() %>% mutate(num_unique=n())
# 
# df_parent2_rep <- df_parent2 %>% filter(nr > num_unique)
# setdiff(df_parent2_rep$cl.x, same_parent_df$cl.x)
# df_parent %>% filter(cl.x=="cc-11")
# 
# # ok now write it out
# 
# 
# df_parent %>% semi_join(df_parent2 %>% filter(nr > num_unique), by=c("cl.x"))
# 
# # do any share accessions?
# mm_counts <- t1_w_info %>% semi_join(mult_acc, by="cl.x")  %>% 
#   select(cl.x, accession) %>%
#   separate_rows(accession, sep="\\|") %>%
#   filter(accession!="list()") %>%
#   group_by(cl.x) %>% mutate(num_accession=n()) %>% 
#   select(cl.x, accession, num_accession) %>% unique() %>% mutate(num_unique=n())
# mm_counts %>% filter(num_accession > num_unique) # none

# synonyms
syn_df <- rbind(t4.2 %>% select(cl, acc),
      t4.2 %>% select(cl, synonym) %>% rename(acc=synonym)) %>%
  arrange(cl, acc) %>%
  unique() %>%
  inner_join(cell_info_df %>% select(primary_accession, derivs, origins, accession), 
             by=c("acc"="primary_accession")) %>%
  separate_rows(accession, sep="\\|") %>%
  filter(accession != "list()") %>%
  separate_rows(derivs, origins, sep="\\|") %>%
  group_by(cl, acc) %>%
  summarize(addl_acc=paste(setdiff(unique(unlist(c(derivs, origins, accession))), c("")), collapse="|"))
t_df <- t_all %>% inner_join(cell_info_df %>% select(primary_accession, derivs, origins, accession), 
                     by=c("acc"="primary_accession")) %>%
  separate_rows(accession, sep="\\|") %>%
  filter(accession != "list()") %>%
  separate_rows(derivs, origins, sep="\\|") %>%
  group_by(cl, acc) %>%
  summarize(addl_acc=paste(setdiff(unique(unlist(c(derivs, origins, accession))), c("")), collapse="|"))

syn_df %>% 
  ungroup() %>%
  separate_rows(addl_acc, sep="\\|") %>%
  filter(addl_acc != "") %>%
  group_by(cl) %>%
  mutate(nr=n()) %>%
  select(cl,  addl_acc, nr) %>%
  unique() %>%
  mutate(num_unique=n()) %>%
  filter(nr > num_unique)


# if two cell lines have the same name and the same parent, then they are mapped to the parent
# if two cell lines have the same name and are different, then they are mapped to both cell lines

# output:
# cl, list_acc, list_dedup_acc, list_syn_acc, list_syn_acc_dedup
t_df2 <- t_df %>% 
  group_by(cl) %>% 
  mutate(overlap=paste(unlist(strsplit(addl_acc, "\\|"))[duplicated(unlist(strsplit(addl_acc, "\\|")))], 
                       collapse="|" )) 
syn_df2 <- syn_df %>%
  ungroup() %>%
  group_by(cl) %>%
  mutate(overlap=paste(unlist(strsplit(addl_acc, "\\|"))[duplicated(unlist(strsplit(addl_acc, "\\|")))], 
                       collapse="|" )) 
  
t_df3 <- t_df2 %>% 
  mutate(res_ids=ifelse(overlap=="", acc, overlap)) 
syn_df3 <- syn_df2 %>%  mutate(res_ids=ifelse(overlap=="", acc, overlap)) 
t_collapse <- t_df3 %>% group_by(cl) %>%
  summarize(list_acc=paste(acc, collapse="|"), list_dedup_acc=paste(unique(res_ids), collapse="|"))
syn_collapse <- syn_df3 %>% group_by(cl) %>% 
  summarize(list_syn=paste(acc, collapse="|"), list_syn_dedup=paste(unique(res_ids), collapse="|"))

syn_w_main <- inner_join(t4.2 %>% select(-synonym) %>% 
                   group_by(cl) %>% summarize(acc=paste(acc, collapse="|")), syn_collapse) 




# multi-mapping data
t_full2 <- t_collapse %>% mutate(list_acc=tolower(list_acc), list_dedup_acc=tolower(list_dedup_acc)) %>%
  group_by(cl) %>% filter(length(strsplit(list_acc, "\\|")[[1]])>1)

# --- add in the extras which have diff punctuation from t3 ---- #
diff_punct <- t3 %>% select(V1) %>% mutate(V1=tolower(V1)) %>% separate(V1, into=c("cl1", "cl2"), sep=", ")
  
  mutate(tolower())
t_full2 %>% head()
diff_punct2 <- diff_punct %>% left_join(t_full2, by=c("cl1"="cl")) %>% select(-cl1) %>% rename(cl=cl2) %>%
  filter(!is.na(list_acc))

multi_map <- rbind(t_full2 %>% ungroup(), diff_punct2 %>% ungroup()) # 408, 62x2 are punctuation

syn_full2 <- syn_w_main %>%  mutate(list_syn=tolower(list_syn), list_syn_dedup=tolower(list_syn_dedup)) %>%
  group_by(cl) %>% filter(length(strsplit(list_syn, "\\|")[[1]])>1)

multi_map %>% write_csv("data/00_db_data/cl_multi_mapping.csv")
syn_full2 %>% write_csv("data/00_db_data/cl_synonym_disambig.csv")

# remove things we're not going to map!
# less than 2 characters, all numeric??
# remove the [ ] part
# ----------------------------- RESTARTING -------------------------- #
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

# -- now need to deal with this! -- #
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
name_parent <- mult_name3%>% filter(overlap!="" & acc != addl_acc) 
name_parent_df <- mult_name4 %>% semi_join(name_parent) # 77
syn_parent_df <- mult_syn4 %>% semi_join(syn_parent) # 739

mult_name4 %>% write_csv("data/00_db_data/cell_line/name_to_acc.csv")
mult_syn4 %>% write_csv("data/00_db_data/cell_line/syn_to_acc.csv")
name_parent_df %>% write_csv("data/00_db_data/cell_line/name_parent.csv")
syn_parent_df %>% write_csv("data/00_db_data/cell_line/syn_parent.csv")


# remove numeric ones



# what are the repeated names?
# length(setdiff(multi_map$cl,cell_lab_name$cl)) # 243 -- I am confused abt this :(
cell_lab_name2 <- cell_lab_name %>% anti_join(multi_map, by=c("cl"))
cell_lab_name3 <- rbind(cell_lab_name2, multi_map %>% select(cl, list_acc) %>% rename(accession=list_acc))

# remove synonyms that are in the list of names
cell_lab_syn2 <- cell_lab_syn %>% anti_join(cell_lab_name3, by=c("synonyms"="cl"))

cell_df <- rbind(dplyr::rename(cell_lab_syn2, cl=synonyms), cell_lab_name3)
cell_df$nwords <- sapply(cell_df$cl, function(x) length(strsplit(x, " ")[[1]]))
cell_df$numchar <- sapply(cell_df$cl, nchar)

#%>%
  #group_by(cl, acc) %>%
  #summarize(list_syn=paste(setdiff(unlist(strsplit(list_acc, "\\|")), unlist(strsplit(acc, "\\|"))), collapse="|"),
  #          list_syn_dedup=paste(setdiff(unlist(strsplit(list_dedup_acc, "\\|")), unlist(strsplit(acc, "\\|"))), collapse="|"))

# 582 cell lines map to a synonym
# remove synonym from


# allele data
alleles <- cell_info_df %>% filter(alleles!= "") # <-- useful for later
alleles %>% select(alleles) %>% unique()
table(alleles$alleles)
# alleles
# 1                 X
# 6               X,Y
# 17            X,Y|X
# 128    Not_detected
# 219           X|X,Y
# 1580 Not_detected|X
# 4578          X,Y|Y
# 5733              Y
