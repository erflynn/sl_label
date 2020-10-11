# Utilities for n-gram mapping
require('tidyverse')
require('rjson')
require('tidytext')
require('tidyverse')
require('fuzzyjoin')
require('stringr')
require('Hmisc')

PUNCT.STR <- "[-|,|\\.|/|#|_|(|)|+|;|:|\t|}|{|\\[|\\]|']"


clean_str <-function(my_str, fill=" ") {
  my_str %>%
    str_replace_all(PUNCT.STR , fill) %>%
    str_squish() %>%
    str_trim() %>%
    tolower()
} 

common_col <- function(dat, col) {
  dat %>% 
    group_by({{col}}) %>% 
    count() %>% 
    arrange(desc(n)) %>%
    ungroup()
}
print_all <- function(dat) print(dat, n=nrow(dat))


# --- code for cleaning metadata --- #
clean_metadata <- function(x, keywords){
  y <- str_replace_all(x, PUNCT.STR, " ");  # substitute punctuation characters
  y <- gsub(keywords, " ", y); # remove specific words
  y <- str_squish(y); # remove repeated whitespace
  z <- str_trim(y) # remove trailing spaces
  return(z)
}

# --- code for labeling ngram data --- #
labelNgram <- function(text_df, ref_df, remove_short=TRUE){
  ref_df2 <- ref_df %>%
    mutate(clean_name=str_trim(str_squish(str_replace_all(name, PUNCT.STR , " ")))) %>%
    unique()
  text_df2 <- text_df %>% 
    mutate(clean_str=str_trim(str_squish(str_replace_all(str, PUNCT.STR , " ")))) %>%
    select(clean_str, str) %>%
    unique()
  
  # separate out the data into unigrams
  sample_unigrams <- text_df2 %>% 
    unnest_tokens(word, clean_str, token=stringr::str_split, pattern = " ") 
  
  if (remove_short){
    sample_unigrams <- sample_unigrams %>% 
      filter(str_detect(word, "[A-z]")) %>% # remove all numeric
      filter(str_length(word) > 3 & is.na(as.numeric(word))) # remove short words or numbers 
  } 
  
  bigrams <- text_df2 %>% 
    unnest_tokens(bigram, clean_str, token = "ngrams",  n = 2)
  
  if (!remove_short){  # remove no+numeric (this is numbering)
    bigrams <- bigrams %>% filter(!str_detect(bigram, "^no [0-9]+"))  
  }
  
  trigrams <- text_df2  %>%
    unnest_tokens(trigram, clean_str, token = "ngrams", n = 3)
  
  # map with a join
  comb_names_uni <- ref_df2 %>% 
    filter(nwords==1) %>%
    inner_join(sample_unigrams, by=c("clean_name"="word")) %>% 
    distinct()
  
  comb_names_bi <- ref_df2 %>% 
    filter(nwords==2) %>%
    inner_join( bigrams, by=c("clean_name"="bigram")) %>% 
    distinct()
  
  comb_names_tri <- ref_df2 %>% 
    filter(nwords==3) %>%
    inner_join(trigrams, by=c("clean_name"="trigram")) %>% 
    distinct()
  
  comb_names <-  bind_rows(comb_names_uni, comb_names_bi, comb_names_tri)
  
  return(comb_names)
}

