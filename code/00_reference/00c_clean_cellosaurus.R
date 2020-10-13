# 
# //TODO - move this cleaning step to later!


require('tidyverse')
require('data.table')


PUNCT.STR <- "[-|,|\\.|/|#|_|(|)|+|;|:|\t]"
PUNCT.STR.SPC <- "[-|,|\\.|/|#|_|(|)|+|;|:|\t| ]"

cell_info_df <- fread("data/00_reference/cellosaurus_df_v2.txt", data.table=FALSE)

cell_dat <- cell_info_df %>% 
  select(cl, primary_accession, synonyms, atcc_acc) %>%
  rename(accession=primary_accession) %>%
  mutate_all(tolower)

cell_lab_name <- cell_dat %>% select(cl, accession) %>% unique()

# add synonyms and atcc data
cell_lab_syn <- cell_dat %>% select(accession, synonyms) %>% 
  separate_rows(synonyms, sep="\\|") %>%
  filter(synonyms != "") %>%
  unique() %>%
  rename(cl=synonyms)

cell_lab_atcc <- cell_dat %>% select(accession, atcc_acc) %>% 
  rename(synonyms=atcc_acc) %>%
  separate_rows(synonyms, sep="\\|") %>%
  filter(synonyms != "") %>%
  unique() %>%
  rename(cl=synonyms)

# remove already present pairs
cell_lab_syn2 <- cell_lab_syn %>% anti_join(cell_lab_name, by=c("accession", "cl")) 
cell_lab_atcc2 <- cell_lab_atcc %>% anti_join(cell_lab_name, by=c("accession", "cl")) 

# remove ones that overlap with names but map to something else
cell_lab_syn3 <- cell_lab_syn2 %>% 
  bind_rows(cell_lab_atcc2) %>% 
  unique() %>%
  anti_join(cell_lab_name, by=c("cl")) # 486

# punctuation overlap
cell_lab_name_punct <- cell_lab_name %>% 
  mutate(cl2=str_replace_all(cl,PUNCT.STR.SPC, ""))  
cell_lab_syn_punct <- cell_lab_syn3 %>% 
  mutate(cl2=str_replace_all(cl,PUNCT.STR.SPC, "")) 

cell_lab_syn_punct2 <- cell_lab_syn_punct %>% 
  anti_join(cell_lab_name_punct, by=c("accession", "cl2"))
overlap_punct <- cell_lab_syn_punct2 %>% semi_join(cell_lab_name_punct, by=c("cl2")) # 319

cell_lab_syn4 <- cell_lab_syn_punct2 %>% anti_join(cell_lab_name_punct, by=c("cl2"))
cell_df <- cell_lab_name_punct %>% select(accession, cl, cl2) %>% bind_rows(cell_lab_syn4)


cell_df2 <- cell_df %>%
  mutate(orig_cl=cl) %>%
  mutate(cl=str_replace_all(cl,PUNCT.STR, "")) %>%
  mutate(cl1=str_replace_all(cl,PUNCT.STR, " ")) %>%
  mutate(cl2=str_replace_all(cl,PUNCT.STR.SPC, "")) %>%
  filter(nchar(cl) >= 3)


discard <- cell_df2$cl[sapply(cell_df2$cl, function(x) 
  str_detect(x, "case|control|extract" ))]
discard_mg  <- cell_df2$cl[sapply(cell_df2$cl, function(x) 
  str_detect(x, "^[0-9| ]+mg$" ))]

other_discard <- c("spindle", "werner", "plate", "dcis", "endo", "smith", "kawasaki", "smoky", 
                   "lima", "star", "la bel", "carb", "scar")
cell_df3 <- cell_df2 %>% filter(!cl %in% c(discard, discard_mg, other_discard))
                     
cell_df3 %>% write_csv("data/00_reference/cell_line/cell_syn_punct_df.csv")
