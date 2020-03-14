require('rjson')
require('tidytext')
require('tidyverse')
require('fuzzyjoin')
require('stringr')
require('Hmisc')
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
sample_metadata <- read_csv(sprintf("data/01_metadata/%s_metadata.csv", prefix))
study_metadata <- read_csv(sprintf("data/01_metadata/%s_experiment_metadata.csv", prefix))

refine_cl_annot <- sample_metadata %>% select(acc, cl_line, part, title) 
  
# clean up the input
text_df <- refine_cl_annot %>% rename(gsm=acc, str=cl_line)
text_df$str <- sapply(text_df$str, function(x){
  y <- str_replace_all(x, "[;\t|(|)]", " ");  ## fix ;\t
  y <- gsub( "\\s+", " ", y);
  z <- gsub(' [[:punct:]]|[[:punct:]] ', ' ', y);  ##remove trailing punctuation
  
  clean_str <- gsub('cells', 'cell', z); # convert cells --> cell
  return(clean_str)
})

cell_df <- read_csv("data/00_db_data/cell_syn_df.csv")

# ---- unigrams ---- #
gsm_unigrams <- text_df %>% 
  unnest_tokens(word, str, token=stringr::str_split, pattern = " ") 
gsm_unigrams2 <- filter(gsm_unigrams, str_length(word) > 3 & is.na(as.numeric(word))) # remove short words or numbers 
gsm_unigrams_char3 <- filter(gsm_unigrams, str_length(word)==3 & is.na(as.numeric(word)))

comb_names <- inner_join(filter(cell_df, nwords==1), gsm_unigrams2, by=c("cl"="word")) %>% distinct()
comb_names_n3 <- inner_join(filter(cell_df, nwords==1), gsm_unigrams_char3, by=c("cl"="word")) %>% distinct()
# these look pretty horrendous! need to pair with outher data

# filter for stopwords
cells <- unique(comb_names$cl )
gsm_map_counts <- gsm_unigrams %>% 
  filter(word %in% cells) %>% 
  group_by(word) %>% summarise(total=n())  %>% 
  arrange(desc(total))

#gsm_map_counts3 <- gsm_unigrams %>% 
#  filter(word %in% comb_names_n3$cl) %>% 
#  group_by(word) %>% summarise(total=n())  %>% 
#  arrange(desc(total))

# selected by looking at gsm_map_counts
cl_stopwords <- c("cancer", "time", "center", "rare", "focus", "peak", "sage", "bona", 
                  "mast", "fisher", "bones", "patches", "madison", "ears", "chance", "cost" )
comb_names2 <- comb_names %>% filter(!cl %in% cl_stopwords)
length(unique(comb_names2$gsm)) # 11251 out of 45036


#comb_data <- right_join(comb_names2[,c("accession", "cl", "gsm")], gsm_data[,c("gsm", "cell_line")])
#write.table(comb_names2[,c("accession", "cl", "gsm")], file="data/tmp_one_word_cell_mapping.txt", row.names=FALSE, sep="\t")


# --- bigrams ---- # 
bigrams <- text_df[,c("gsm", "str")] %>% unnest_tokens(bigram, str, token = "ngrams",  n = 2)
bigrams_separated <- bigrams %>% separate(bigram, c("word1", "word2"), sep = " ")
expanded_stop_words <- c(stop_words$word, "human", "hospital")
bigrams_filt <- bigrams_separated %>% 
  filter(is.na(as.numeric(word1)))  %>% filter(is.na(as.numeric(word2))) %>%
  filter(!word1 %in% expanded_stop_words) %>%
  filter(!word2 %in% expanded_stop_words)
# filter for stopwords, numeric, etc
bigram_counts <- bigrams_filt %>% 
  count(word1, word2, sort = TRUE)

bigrams_united <- bigrams_filt %>%
  unite(bigram, word1, word2, sep = " ")

#cell_two_words <- rbind(filter(cell_df2, nwords==2), cell_df2_short)
cells2 <- cell_df %>% filter(nwords==2)
bigram_map_counts <- bigrams_united %>% filter(bigram %in% cells2) %>% group_by(bigram) %>% summarise(total=n())  %>% arrange(desc(total))
mapped_two <- inner_join(bigrams_united, cells2, c("bigram"="cl"))


# --- trigrams --- #   
trigrams <- text_df[,c("gsm", "str")] %>% unnest_tokens(trigram, str, token = "ngrams", n = 3)
#cell_more <- filter(cell_df, nwords>4) %>% unnest_tokens(cl, cl, token="ngrams", n=3)
cell_three_words <- filter(cell_df, nwords==3)
mapped_three <- inner_join(trigrams, cell_three_words, c("trigram"="cl"))
trigram_counts <- mapped_three %>% group_by(trigram) %>% summarise(total=n())  %>% arrange(desc(total))


# --- what to do about bigger?? 3 words --- #


# ---- which contain ATCCC ? map those ---- #
# atcc
text_df %>% filter(stringr::str_detect(str, "atcc"))


# ---- deal with dashes/spaces ---- #

# partial match
cell_no_dash <- cell_df %>% mutate(cl=str_replace_all(cl, "-", " "))
cell_no_spc <- cell_df %>% mutate(cl=str_replace_all(cl, " ", ""))
cell_no_dash_spc <- cell_df %>% mutate(cl=str_replace_all(cl, "[ |-]", ""))

# what about ones that multi-map??

# --- look at 100 that map and 100 that don't --- #
no_map <- sample_metadata %>% 
  anti_join(comb_names , by=c("acc"="gsm")) %>% 
  filter(!is.na(cl_line)) %>%
  group_by(cl_line) %>% 
  count() %>%
  arrange(desc(n))

no_map %>% filter(str_detect(cl_line, "cell"))

### put together all the mappings
mapped_dat <- do.call(rbind, list(comb_names2[,c("accession", "cl", "gsm")],
                                  mapped_two %>% select(gsm, bigram, accession) %>% dplyr::rename(cl=bigram),
                                  mapped_three %>% select(gsm, trigram, accession) %>% dplyr::rename(cl=trigram)))

# write out
mapped_all <- mapped_dat %>% group_by(gsm) %>% 
  mutate( cl=paste(unique(strsplit(cl, ",")), collapse=";"), 
          accession=paste(unique(strsplit(accession, ",")), collapse=";")) 
mapped_all2 <- mapped_all %>% ungroup() %>% unique()
mapped_all2 %>% write_csv(sprintf("data/%s_acc_to_cl.csv", prefix))
