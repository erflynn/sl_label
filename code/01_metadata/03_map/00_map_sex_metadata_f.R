
require('tidyverse')
PUNCT.STR <- "[-|,|\\.|/|#|_|(|)|+|;|:|\t|}|{|\\[|\\]|']"


# --- manually selected lists of keywords --- #
hc_missing_keywords <- c("missing", "unknown", "undetermined")
missing_keywords <- c( "missing", "unknown", "na", "n a", "not available", 
                       "",  "undetermined", "other", "u",
                       "not specified", "not provided", "not collected", "nk", 
                       "n k", "n d", "nd", "no data", "not determined", 
                       "not applicable", "none", "?", "indeterminate",
                       "undetermine",
                       "unkown")
mixed_keywords <- c("mixed", "pooled", "mix", "mixture", "both")
f_keywords <- c("female", "females", "xx", "fem", "f", "woman", "women", 
                "girl", "femal", 
                "femile", "fimale", "femlale")
m_keywords <- c("male", "males", "xy", "mal", "m", "man", "men", "boy")
hc_m_keywords <- c("male", "males", "man", "men")
hc_f_keywords <- c("female", "females", "woman", "women")

  
clean_metadata <- function(x, keywords){
  y <- str_replace_all(x, PUNCT.STR, " ");  # substitute punctuation characters
  y <- gsub(keywords, " ", y); # remove the word sex
  y <- str_squish(y); # remove repeated whitespace
  z <- str_trim(y) # remove trailing spaces
  return(z)
}

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

