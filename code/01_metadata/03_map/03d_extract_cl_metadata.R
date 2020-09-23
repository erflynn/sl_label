# code for exctracting cell line data
require('data.table')
require('tidyverse')
source("code/01_metadata/03_map/00_mapping_utils.R")

clean_str <-function(my_str, fill=" ") {
  my_str %>%
    str_replace_all(PUNCT.STR , fill) %>%
    str_squish() %>%
    str_trim() %>%
    tolower()
} 

load("data/01_metadata/all_sample_attrib_clean.RData") # --> all_attrib_clean


# --- divide into human + mouse data --- #
cell_info_df <- fread("../labeling/geo2drug/data/00_db_data/cellosaurus_df_v2.txt", data.table=FALSE)
cell_dat <- cell_info_df %>% select(species, primary_accession, cl) %>%
  mutate(across(c(primary_accession, cl), tolower)) %>%
  rename(accession=primary_accession, cl_name=cl) %>%
  group_split(species) 
human_cl <- cell_dat[[1]]
mouse_cl <- cell_dat[[2]]


# --- deal with the cell line data --- #
cell_df2 <- read_csv("data/00_db_data/cell_line/cell_syn_punct_df.csv")
cell_df2 <- cell_df2 %>% 
  mutate(cl3=clean_str(orig_cl, fill="")) # add other cleaning
cell_df_nodash <- do.call(rbind,
                          list(cell_df2 %>% select(accession, cl3) %>% rename(cl=cl3), 
                               cell_df2 %>% select(accession, cl2) %>% rename(cl=cl2) ,
                               cell_df2 %>% select(accession, cl1) %>% rename(cl=cl1),
                               cell_df2 %>% select(accession, cl))) %>% unique()


cl_stopwords <- c("cancer", "time", "center", "rare", "focus", "peak", "sage", "bona", 
                  "mast", "fisher", "bones", "patches", "madison", "ears", 
                  "chance", "cost", "wicell", "time", "huvec", "human lung", "129s6")

cell_df3 <- cell_df_nodash %>% 
  filter(! cl %in% cl_stopwords) %>%
  mutate(
    nwords=length(strsplit(cl, " ")[[1]]),
    numchar=nchar(cl)) 


# separate out
human_cl2 <- inner_join(human_cl, cell_df3, by="accession") %>%
  select(-species) %>% unique()
mouse_cl2 <- inner_join(mouse_cl, cell_df3, by="accession") %>%
  select(-species) %>% unique()

# //TODO: age: p100, p162, etc is bad
mouse_cl2 %>% 
  filter(!str_detect(cl, "^p[0-9]+$")) # probably remove these

mouse_strains <- read_csv("data/00_db_data/jackson_lab_strains.csv", col_types="c")
mouse_strains2 <- mouse_strains %>% mutate(Strain=tolower(Strain))
mouse_strains3 <- mouse_strains2 %>% 
  mutate(Strain1=clean_str(Strain, fill=" ")) %>%
  mutate(Strain2=clean_str(Strain, fill="")) %>%
  group_by(Strain) %>%
  mutate(Strain4=str_split(Strain1, " ")[[1]][[1]]) %>%
  ungroup() 
intersect(mouse_strains3$Strain4, mouse_cl2$cl)

# --- separate out the other data --- #
comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", 
                          col_types="cccccccdld")

human_s <- comb_metadata %>% filter(organism=="human")
mouse_s <- comb_metadata %>% filter(organism=="mouse")

# --- get the appropriate attribute data --- #
cell_line <- all_attrib_clean %>% 
  filter(str_detect(key_clean, "cell"),
         str_detect(key_clean, "line"))

cell_line2 <- all_attrib_clean %>% 
  anti_join(cell_line) %>%
  filter(str_detect(value_clean, "\\bcell"),
         str_detect(value_clean, "\\bline"))

# broader
cell_data <- all_attrib_clean %>% 
  filter((str_detect(key_clean, "\\bcell") |
            str_detect(value_clean, "\\bcell")) & 
           !str_detect(key_clean, "tissue")) 
# mess up -- key: tissue, value:lncap



summarizeMapped <- function(df, label){
  df %>% 
    distinct(sample_acc, accession) %>%
    group_by(sample_acc) %>%
    summarise(accession=paste(accession, collapse=";")) %>%
    ungroup() %>%
    mutate(type=label)
}


mapClDat <- function(list_samples, ref_cl){
  h_cell_line <- cell_line %>% 
    semi_join(list_samples, by="sample_acc")  %>%
    mutate(value=tolower(value)) %>%
    filter(is.na(as.numeric(value_clean))) 
  h_cell_line2 <- cell_line2 %>% 
    semi_join(list_samples, by="sample_acc") %>%
    mutate(value=tolower(value))
  h_cell_data <- cell_data %>% 
    semi_join(list_samples, by="sample_acc") 
  
  # ---- exact match ----- #
  ref_cl2 <- ref_cl %>% 
    distinct(accession, cl_name) %>%
    group_by(cl_name) %>%
    mutate(nwords=length(str_split(cl_name, " ")[[1]]), 
           numchar=nchar(cl_name)) %>%
    ungroup()
  num_samples <- h_cell_line %>% 
    nrow() # human: 99,426
  num_cl_fields <- h_cell_line %>% 
    filter(is.na(as.numeric(value_clean))) %>% 
    distinct(value_clean) %>% 
    nrow() # human: 7,587
  print(num_samples)
  print(num_cl_fields)
  
  
  exact_match <- h_cell_line %>% 
    distinct(value) %>% 
    filter(is.na(as.numeric(value))) %>%
    inner_join(ref_cl2, by=c("value"="cl_name")) 
  # TODO - some of the duplicate mappings are case-sensitive :(  
  
  exact_match2 <- h_cell_line %>% 
    distinct(value_clean) %>% 
    inner_join(ref_cl2, by=c("value_clean"="cl_name")) 
  
  exact_match3 <- h_cell_line %>% 
    distinct(value) %>% 
    filter(is.na(as.numeric(value))) %>%
    inner_join(ref_cl, by=c("value"="cl")) 
  # TODO check this
  
  exact_match4 <- h_cell_line %>%
    mutate(value_clean2=clean_str(value, fill="")) %>%
    distinct(value, value_clean2) %>%
    filter(is.na(as.numeric(value_clean2))) %>%
    inner_join(ref_cl, by=c("value_clean2"="cl")) 
  
  exact_mapped <- h_cell_line %>% 
    inner_join(exact_match, by="value")
  exact_mapped2 <- h_cell_line %>% 
    anti_join(exact_match, by="value") %>%
    anti_join(exact_mapped, by="sample_acc") %>%
    inner_join(exact_match2, by="value_clean")
  
  exact_mapped3 <- h_cell_line %>% 
    anti_join(exact_mapped, by="sample_acc") %>%
    anti_join(exact_mapped2, by="sample_acc") %>%
    inner_join(exact_match3, by="value")
  exact_mapped4 <- h_cell_line %>% 
    anti_join(exact_mapped, by="sample_acc") %>%
    anti_join(exact_mapped2, by="sample_acc") %>%
    anti_join(exact_mapped3, by="sample_acc") %>%
    inner_join(exact_match4, by="value")
  
  mapped_dat <- exact_mapped %>% summarizeMapped("exact") %>%
    bind_rows(exact_mapped2 %>% summarizeMapped("exact2"),
              exact_mapped3 %>% summarizeMapped("exact3"),
              exact_mapped4 %>% summarizeMapped("exact4")) 
  
  
  # --- map using fields --- #
  h_cell_line_filt <- h_cell_line %>% 
    anti_join(exact_match, by="value") %>%
    anti_join(exact_match2, by="value_clean") %>%
    anti_join(mapped_dat, by="sample_acc") 
  
  h_map1 <- labelNgram(h_cell_line_filt %>% 
                         rename(str=value_clean, 
                                sample=sample_acc) %>% 
                         mutate(orig_str=value) %>% 
                         select(sample, orig_str, str), 
                       ref_cl %>% rename(name=cl))
  
  h_cell_line_mapped <- h_cell_line_filt %>%
    inner_join(h_map1, by=c("value_clean"="str"))  
    
  mapped_dat_f <- h_cell_line_mapped %>% summarizeMapped("cell_line_field")
  mapped_dat2 <- mapped_dat %>% bind_rows(mapped_dat_f)
  num_mapped <- mapped_dat2 %>% distinct(sample_acc) %>% nrow() # 74,140
  num_str_mapped <- h_cell_line %>% 
    semi_join(mapped_dat2, by="sample_acc") %>%
    distinct(value_clean) %>%
    nrow()# 4,884
  print(num_mapped)
  print(num_str_mapped)
  h_cell_line2_filt <- h_cell_line2 %>% anti_join(mapped_dat2, by="sample_acc")
  
  h_map2 <- labelNgram(h_cell_line2_filt %>% 
                         rename(str=value_clean, sample=sample_acc) %>% 
                         mutate(orig_str=value) %>% 
                         select(sample, orig_str, str), 
                       ref_cl %>% rename(name=cl))
  h_map2 %>% distinct(accession) %>% nrow() 
  
  h_cell_line_mapped2 <- h_cell_line2_filt %>% 
    inner_join(h_map2, by=c("value_clean"="str")) 
  
  mapped_dat_f2 <- h_cell_line_mapped2 %>% summarizeMapped("cell_line_mention")
  
  mapped_dat3 <- mapped_dat2 %>% bind_rows(mapped_dat_f2)
  cell_data_f <- h_cell_data %>% anti_join(mapped_dat3, by="sample_acc")
  
  # --- weaker mapping --- #
  map3 <- labelNgram(cell_data_f %>% 
                       rename(str=value_clean, sample=sample_acc) %>% 
                       mutate(orig_str=value) %>% 
                       select(sample, orig_str, str), 
                     ref_cl %>% rename(name=cl))

  h_cell_line_mapped3 <- cell_data_f %>% 
    inner_join(map3, by=c("value_clean"="str"))  
  mapped_dat_f3 <- h_cell_line_mapped3 %>% summarizeMapped("other")
  mapped_dat4 <- mapped_dat3 %>% bind_rows(mapped_dat_f3)
  
  return(mapped_dat4)
}
mapped_human <- mapClDat(human_s, human_cl2)
# 99426, 7587 --> 74140, 4884
mapped_mouse <- mapClDat(mouse_s, mouse_cl2)
# 17061, 1186 --> 8433, 309

table(mapped_human$type)
table(mapped_mouse$type)

# -- compare to MetaSRA -- #
require('rjson')
cvcl_map <- fromJSON(file="../metasra_cvcl_mappings.json") 

cvcl_df <- do.call(rbind, lapply(cvcl_map, function(x) 
  paste(x$mapped_terms, collapse=";")))

cell_mapped <- mapped_df %>% filter(str_detect(term_id, "CVCL")) %>%
  mutate(term_id2=str_replace_all(tolower(term_id), ":", "_"))
cl_mapped2 <- cell_mapped %>% left_join(samp_to_runs, by=c("sample_accession"="sample")) %>%
  dplyr::rename(acc=run) 
head(cl_mapped2)

human_rnaseq <- comb_metadata %>% filter(data_type=="rnaseq", organism=="human", num_reads > 100000, present)
mapped_human_rnaseq <- mapped_human %>% semi_join(human_rnaseq)
cl_mapped3 <- cl_mapped2 %>% semi_join(human_rnaseq, by=c("acc"="sample_acc"))
cl_mapped3 %>% anti_join(mapped_human_rnaseq, by=c("acc"="sample_acc")) %>%
  sample_n(10)
length(setdiff(cl_mapped3$acc, mapped_human_rnaseq$sample_acc)) # 6611
length(setdiff(mapped_human_rnaseq$sample_acc, cl_mapped3$acc)) # 15459

both_cl <- cl_mapped2 %>% inner_join(mapped_human_rnaseq, by=c("acc"="sample_acc"))
both_cl %>% filter(term_id2==accession) %>% nrow() # 14390/18821 (76.5%)
multi_map <- both_cl %>% filter(term_id2!=accession & str_detect(accession, ";")) 
some_overlap <- multi_map %>% filter(str_detect(accession,term_id2)) # 629
no_overlap <- multi_map %>% filter(!str_detect(accession,term_id2)) # 326

bind_rows(mapped_human, mapped_mouse) %>% 
  write_csv("data/02_labeled_data/cell_line_mapping.csv")

# metrics:
#  [x] compare to previous
#  [x] fraction of cell line fields mapping
#  [x] compare to metaSRA
#  - 100 randomly selected that do not map but contain a cell word 
#  - 100 randomly selected that do map