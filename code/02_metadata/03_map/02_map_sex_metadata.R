# code for mapping sample sex metadata

require('tidyverse')
source("code/01_metadata/03_map/00_mapping_utils.R")

# --- manually selected lists of keywords --- #
hc_missing_keywords <- c("missing", "unknown", "undetermined")
missing_keywords <- c( "missing", "unknown", "na", "n a", "not available", 
                       "",  "undetermined", "other", "u",
                       "not specified", "not provided", "not collected", "nk", 
                       "n k", "n d", "nd", "no data", "not determined", 
                       "not applicable", "none", "?", "indeterminate",
                       "undetermine", "unkown")
mixed_keywords <- c("mixed", "pooled", "mix", "mixture", "both")
f_keywords <- c("female", "females", "xx", "fem", "f", "woman", "women", 
                "girl", "femal", 
                "femile", "fimale", "femlale")
m_keywords <- c("male", "males", "xy", "mal", "m", "man", "men", "boy")
hc_m_keywords <- c("male", "males", "man", "men")
hc_f_keywords <- c("female", "females", "woman", "women")



# --- function for mapping sex labels to one of f, m, mixed, unknown --- #
map_sex_metadata <- function(dat, col){
  dat_unique <- dat %>% distinct({{col}}) %>% rename(sex={{col}})
  dat_unique <- dat_unique %>% 
    group_by(sex) %>%
    mutate(reform_sex=clean_metadata(sex,"sex|gender"))
  
  # convert each label to one of these
  dat_unique <- dat_unique %>% mutate(
    annot_sex=case_when(
      is.na(reform_sex) ~ "unknown",
      reform_sex %in% missing_keywords ~ "unknown",
      reform_sex %in% mixed_keywords ~ "mixed",
      reform_sex %in% f_keywords ~ "female",
      reform_sex %in% m_keywords ~ "male",
      reform_sex=="1" ~ "male", # these are only exact matches
      reform_sex=="0" ~ "female", # only exact matches
      str_detect(reform_sex, "^[f]+[^m]$") ~ "female",
      str_detect(reform_sex, "^[m]+[^f]$") ~ "male",
      str_detect(reform_sex, "^[fm]+$") & 
        str_detect(reform_sex, "f") & 
        str_detect(reform_sex,"m") ~ "mixed",
      TRUE ~ ""
    )
  )
  
  # deal with multiple words
  dat_unique2 <- dat_unique %>% 
    mutate(tokens=stringr::str_split(reform_sex, pattern=" "))
  
  dat_unique3 <- dat_unique2 %>%
    mutate(annot_sex2=case_when(    
      #(length(intersect(tokens[[1]], mixed_keywords)) > 0) ~ "mixed",
      (length(intersect(tokens[[1]], hc_f_keywords)) > 0  &
         length(intersect(tokens[[1]], hc_m_keywords)) > 0) ~ "mixed",
      (length(intersect(tokens[[1]], hc_f_keywords)) > 0)  ~ "female",
      (length(intersect(tokens[[1]], hc_m_keywords)) > 0)  ~ "male",
      (length(intersect(tokens[[1]], hc_missing_keywords)) > 0)  ~ "unknown",
      TRUE ~ ""
    ))
  
  # collapse the annotated sex labels
  dat_unique4 <- dat_unique3 %>% 
    mutate(mapped_sex=case_when(
      annot_sex != "" ~ annot_sex,
      annot_sex2 != "" ~ annot_sex2,
      TRUE ~ ""
    )) %>%
    select(sex, mapped_sex) %>%
    ungroup()
  
  return(dat_unique4)
}


# load the data
load("data/01_metadata/all_sample_attrib_clean.RData") # --> all_attrib_clean


#  grab attribute pairs that are relevant
sg_sample_attr <- all_attrib_clean %>% 
  filter(str_detect(key_clean, "\\bsex\\b|\\bgender\\b"))
not_sg <- all_attrib_clean %>% 
  anti_join(sg_sample_attr, by="sample_acc") 
rescue_sg <- not_sg  %>% 
  filter(!str_detect(key_clean, "strain|genetic|recipient|donor|child|sex chromosome complement")) %>% 
  filter(str_detect(value_clean, "male"))

# sample, key, value, type_key (std/rescue), normalized_value
sg_mapped <- sg_sample_attr %>% map_sex_metadata(value_clean)
sg_mapped2 <- rescue_sg %>% map_sex_metadata(value_clean)
sg_tab <- sg_sample_attr %>% 
  mutate(type_key="std") %>%
  bind_rows(rescue_sg %>% mutate(type_key="rescue")) 
sg_tab2 <- sg_tab %>% 
  left_join(sg_tab %>% map_sex_metadata(value_clean), 
            by=c("value_clean"="sex"))

# de-duplicate and condense such that every row is a sample
sg_tab3 <- sg_tab2 %>% unique() %>% filter(!str_detect(key_clean, "recipient|donor|child|sex chromosome complement")) 

mult_row <- sg_tab3 %>%  group_by(sample_acc) %>% count() %>% arrange(desc(n)) %>% filter(n > 1)
single_row <- sg_tab3 %>% anti_join(mult_row) %>% select(sample_acc, key_clean, value_clean, type_key, mapped_sex)

# for multi-row, combine and label with the non-unknown label
# if maps to multiple sex labels, label as mixed
mult_row_condensed <- sg_tab3 %>% semi_join(mult_row) %>%
  filter(mapped_sex != "unknown") %>%
  group_by(sample_acc) %>%
  summarise(key_clean=paste(unique(key_clean), collapse=";"),
            value_clean=paste(unique(value_clean), collapse=";"),
            type_key=paste(unique(type_key), collapse=";"),
            mapped_sex=paste(unique(sort(mapped_sex)), collapse=";"))  %>%
  mutate(mapped_sex=ifelse(mapped_sex=="female;male", "mixed", mapped_sex))

sg_tab_mapped <- single_row %>% bind_rows(mult_row_condensed)
stopifnot(nrow(sg_tab_mapped)==length(sg_tab_mapped$sample_acc))

sg_tab_mapped %>% write_csv("data/01_metadata/mapped_sl_all.csv")

