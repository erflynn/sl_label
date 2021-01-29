# clean_sex_metadata.R
# E Flynn
#
# Code to grab sex annotations that don't fall into "male" or "female"
# These are then mapped using manually constructed rules to male, female, mixed, or unknown.


require('tidyverse')
require('tidytext')

options(stringsAsFactors = FALSE)

# microarray
human_metadata <- read.csv("data/01_metadata/human_metadata.csv")
mouse_metadata <- read.csv("data/01_metadata/mouse_metadata.csv")
rat_metadata <- read.csv("data/01_metadata/rat_metadata.csv")

# rna-seq
human_rnaseq_metadata <- read.csv("data/01_metadata/human_rnaseq_sample_metadata.csv")
mouse_rnaseq_metadata <- read.csv("data/01_metadata/mouse_rnaseq_sample_metadata.csv")
rat_rnaseq_metadata <- read.csv("data/01_metadata/rat_rnaseq_sample_metadata.csv")

# grab the possible sex labels
list_metadata <- list(human_metadata, mouse_metadata, rat_metadata,
                      human_rnaseq_metadata, mouse_rnaseq_metadata, rat_rnaseq_metadata)
list_sex_lab <- lapply(list_metadata, function(x) x %>% select(sex) %>% unique())

# manually selected lists of keywords
hc_missing_keywords <- c("missing", "unknown", "undetermined")
missing_keywords <- c( "missing", "unknown", "na", "n a", "not available", 
                       "-", "--", "-", "", ".", "undetermined", "other",
                       "not specified", " ", "not provided", "not collected", "nk", 
                       "n k", "n d", "nd", "no data", "not determined", 
                       "not applicable", "none", "?", "indeterminate",
                       "unkown")
mixed_keywords <- c("mixed", "pooled", "mix", "mixture", "both")
f_keywords <- c("female", "females", "xx", "fem", "f", "woman", "women", 
                "girl", "femal", 
                "femile", "fimale", "femlale")
m_keywords <- c("male", "males", "xy", "mal", "m", "man", "men", "boy")

unique_lab <- do.call(rbind, list_sex_lab)

# clean punctuation
unique_lab$reform_sex <- lapply(unique_lab$sex, function(x){
  y <- str_replace_all(x, "[;|\\(|\\)|\\/|\\_|,|:|\\.|\\+|\\-|\\']", " ");  # substitute punctuation characters
  y <- gsub("sex", " ", y); # remove the word sex
  y <- str_squish(y); # remove repeated whitespace
  z <- str_trim(y) # remove trailing spaces
  return(z)
})

# convert each label to one of these
unique_lab2 <- unique_lab %>% mutate(
  annot_sex=case_when(
    is.na(reform_sex) ~ "unknown",
    reform_sex %in% missing_keywords ~ "unknown",
    reform_sex %in% mixed_keywords ~ "mixed",
    reform_sex %in% f_keywords ~ "female",
    reform_sex %in% m_keywords ~ "male",
    reform_sex=="1" ~ "male", # these are only exact matches
    reform_sex=="0" ~ "female", # only exact matches
    TRUE ~ ""
  )
)

# deal with multiple words
unique_lab2$tokens <- sapply(unique_lab2$reform_sex, stringr::str_split, pattern=" ")

unique_lab2$annot_sex2 <- 
  sapply(unique_lab2$tokens , function(tokens){
  case_when(    
    (length(intersect(tokens, mixed_keywords)) > 0) ~ "mixed",
    (length(intersect(tokens, f_keywords)) > 0  &
      length(intersect(tokens, m_keywords)) > 0) ~ "mixed",
    (length(intersect(tokens, f_keywords)) > 0)  ~ "female",
    (length(intersect(tokens, m_keywords)) > 0)  ~ "male",
    (length(intersect(tokens, hc_missing_keywords)) > 0)  ~ "unknown",
    TRUE ~ ""
)})
# unique_lab2 %>% filter(annot_sex2=="" & annot_sex=="") %>% View()
# unique_lab2 %>% filter(annot_sex2=="mixed" | annot_sex=="mixed") %>% View()
# unique_lab2 %>% filter(annot_sex2=="male" | annot_sex=="male") %>% View()
# unique_lab2 %>% filter(annot_sex2=="unknown" | annot_sex=="unknown") %>% View()

# collapse the annotated sex labels
unique_lab3 <- unique_lab2 %>% mutate(
  mapped_sex=case_when(
    annot_sex != "" ~ annot_sex,
    annot_sex2 != "" ~ annot_sex2,
    TRUE ~ "unknown")
  )
unique_lab4 <- unique_lab3 %>%
  select(sex, mapped_sex) %>%
  unique()

# write out the mappings that have been made  
unique_lab4 %>% write_csv("data/alt_sex_annot.csv")

# add in the mapped sex labels
list_metadata <- list(human_metadata, mouse_metadata, rat_metadata,
                      human_rnaseq_metadata, mouse_rnaseq_metadata, rat_rnaseq_metadata)
list_metadata2 <- lapply(list_metadata, function(x) 
  x %>% select(acc, sex) %>% left_join(unique_lab4) )

# now write them all out
list_metadata2[[1]] %>% unique() %>% write_tsv("data/01_metadata/human_microarray_metadata_sex.tsv")
list_metadata2[[2]] %>% unique() %>% write_tsv("data/01_metadata/mouse_microarray_metadata_sex.tsv")
list_metadata2[[3]] %>% unique() %>% write_tsv("data/01_metadata/rat_microarray_metadata_sex.tsv")
list_metadata2[[4]] %>% unique() %>% write_tsv("data/01_metadata/human_rnaseq_metadata_sex.tsv")
list_metadata2[[5]] %>% unique() %>% write_tsv("data/01_metadata/mouse_rnaseq_metadata_sex.tsv")
list_metadata2[[6]] %>% unique() %>% write_tsv("data/01_metadata/rat_rnaseq_metadata_sex.tsv")
