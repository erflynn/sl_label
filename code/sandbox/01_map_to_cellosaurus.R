require('tidyverse')
require('rjson')
require('tidytext')
require('tidyverse')
require('fuzzyjoin')
require('stringr')
require('Hmisc')

# // TODO use labelNgram

options(stringsAsFactors = FALSE)
# map both study + sample

args <- commandArgs(trailingOnly=TRUE)

prefix <- args[1]
data_type <- args[2]


PUNCT.STR <- "[-|,|\\.|/|#|_|(|)|+|;|:|\t]"
PUNCT.STR.SPC <- "[-|,|\\.|/|#|_|(|)|+|;|:|\t| ]"


mapTextCl <- function(text_df, cell_df, three_l=TRUE){
  cell_df <- cell_df %>% mutate(
    nwords=length(strsplit(cl, " ")[[1]]),
    numchar=nchar(cl)) 
  # ---- unigrams ---- #
  gsm_unigrams <- text_df %>% 
    unnest_tokens(word, str, token=stringr::str_split, pattern = " ") %>%
    filter( is.na(as.numeric(word))) # remove numbers
  gsm_unigrams2 <- filter(gsm_unigrams, 
                          str_length(word) > 3) # remove short words 
  gsm_unigrams_char3 <- filter(gsm_unigrams, 
                               str_length(word)==3)
  
  comb_names <- inner_join(filter(cell_df, nwords==1), 
                           gsm_unigrams2, by=c("cl"="word")) %>% 
    distinct()
  comb_names_n3 <- inner_join(filter(cell_df, nwords==1), 
                              gsm_unigrams_char3, by=c("cl"="word")) %>% 
    distinct()
  # these look pretty horrendous! need to pair with outher data
  
  # filter for stopwords
  cells <- unique(comb_names$cl)
  gsm_map_counts <- gsm_unigrams %>% 
    filter(word %in% cells) %>% 
    group_by(word) %>% summarise(total=n())  %>% 
    arrange(desc(total))
  
  gsm_unigrams_char3 %>% group_by(word) %>% count() %>% arrange(desc(n))
  char3_stop <- c("the", "and", "sum", "age", "hey", "for", "wbs", "mel", "ctr", "lcl")
  comb_names_n3.2 <- comb_names_n3 %>% filter(!cl %in% char3_stop)
  # <--- I don't think we should use.... 
  # or if we do it needs to be 1word? or paired or something?---> #
  
  # selected by looking at gsm_map_counts
  cl_stopwords <- c("cancer", "time", "center", "rare", "focus", "peak", "sage", "bona", 
                    "mast", "fisher", "bones", "patches", "madison", "ears", "chance", "cost" )
  comb_names2 <- comb_names %>% filter(!cl %in% cl_stopwords)
  length(unique(comb_names2$orig_str)) # 11251 out of 45036
  if (three_l){
    comb_names2 <- bind_rows(comb_names2, comb_names_n3.2)
  }
  
  # --- bigrams ---- # 
  bigrams <- text_df[,c("orig_str", "str")] %>% unnest_tokens(bigram, str, token = "ngrams",  n = 2)
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
  mapped_two %>% group_by(bigram) %>% count() %>% arrange(desc(n))
  
  # --- trigrams --- #   
  trigrams <- text_df[,c("orig_str", "str")] %>% unnest_tokens(trigram, str, token = "ngrams", n = 3)
  #cell_more <- filter(cell_df, nwords>4) %>% unnest_tokens(cl, cl, token="ngrams", n=3)
  cell_three_words <- filter(cell_df, nwords==3)
  mapped_three <- inner_join(trigrams, cell_three_words, c("trigram"="cl"))
  trigram_counts <- mapped_three %>% group_by(trigram) %>% 
    summarise(total=n())  %>% arrange(desc(total))
  
  ### put together all the mappings
  mapped_dat <- do.call(rbind, list(comb_names2[,c("accession", "cl", "orig_str")],
                                    mapped_two %>% select(orig_str, bigram, accession) %>% 
                                      dplyr::rename(cl=bigram),
                                    mapped_three %>% select(orig_str, trigram, accession) %>% 
                                      dplyr::rename(cl=trigram)))
  
  # write out
  mapped_all <- mapped_dat %>% group_by(orig_str) %>% 
    mutate( cl=paste(unique(strsplit(cl, ",")), collapse=";"), 
            accession=paste(unique(strsplit(accession, ",")), collapse=";")) 
  mapped_all2 <- mapped_all %>% ungroup() %>% unique()
  return(mapped_all2)
}

# options: 
# - all punctuation converted to ""
# - all punctuation converted to " "
# - all spaces + punctuation converted to ""

if (data_type=="compendia"){
  sample_metadata <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", prefix))
  study_metadata <- read.csv(sprintf("data/01_metadata/%s_experiment_metadata.csv", prefix))
  
} else {
  sample_metadata <- read.csv(sprintf("data/01_metadata/%s_%s_sample_metadata.csv", prefix, data_type))
  study_metadata <- read.csv(sprintf("data/01_metadata/%s_%s_experiment_metadata.csv", prefix, data_type))
  
}

refine_cl_annot <- sample_metadata %>% select(acc, cl_line, part, title)  %>%
  filter(!is.na(cl_line) & cl_line != "" & cl_line != "--")
refine_cl_annot %>% filter(str_detect(cl_line, "primary")) %>% select(cl_line) %>% unique()
refine_cl_annot %>% filter(str_detect(cl_line, "stem")) %>% select(cl_line) %>% unique()
refine_cl_annot %>% filter(str_detect(cl_line, "culture")) %>% select(cl_line) %>% unique()

study_text_df <- study_metadata %>% 
  rename(gse=study_acc) %>% 
  group_by(gse) %>%
  mutate(orig_str=paste(c(title, description), collapse=" ")) %>%
  select(-title, -description) %>%
  mutate(orig_str=tolower(orig_str)) %>%
  mutate(str=str_trim(str_squish(str_replace_all(orig_str, PUNCT.STR , "")))) %>%
  mutate(str1=str_trim(str_squish(str_replace_all(orig_str, PUNCT.STR , " "))))

# clean up the input
text_df <- refine_cl_annot %>% rename(gsm=acc, orig_str=cl_line) %>%
  mutate(orig_str=tolower(orig_str)) %>%
  group_by(orig_str) %>%
  filter(nchar(orig_str) >= 3 & is.na(as.numeric(orig_str))) %>%
  ungroup() %>%
  mutate(str=str_trim(str_squish(str_replace_all(orig_str, PUNCT.STR , "")))) %>%
  mutate(str1=str_trim(str_squish(str_replace_all(orig_str, PUNCT.STR , " ")))) %>%
  mutate(str2=str_trim(str_squish(str_replace_all(orig_str, PUNCT.STR.SPC , "")))) 



text_df_title <- bind_rows(sample_metadata %>% select(acc, part) %>% rename(orig_str=part),
                           sample_metadata %>% select(acc, title) %>% rename(orig_str=title)) %>%
  rename(gsm=acc) %>%
  mutate(orig_str=tolower(orig_str)) %>%
  group_by(orig_str) %>%
  filter(nchar(orig_str) >= 3 & is.na(as.numeric(orig_str))) %>%
  ungroup() %>%
  mutate(str1=str_trim(str_squish(str_replace_all(orig_str, PUNCT.STR , " ")))) 

text_df2 <- text_df %>% 
  select(orig_str, str, str1, str2) %>% 
  unique() %>% 
  filter(!str_detect(orig_str, "c57|bl6|balb")) 

text_df_title2 <- text_df_title %>% 
  select(orig_str, str1) %>% 
  unique() %>% 
  filter(!str_detect(orig_str, "c57|bl6|balb"))

cell_df2 <- read_csv("data/00_db_data/cell_line/cell_syn_punct_df.csv")
cell_df_nodash <- do.call(rbind,
                          list(cell_df2 %>% select(accession, cl2) %>% rename(cl=cl2) ,
                               cell_df2 %>% select(accession, cl1) %>% rename(cl=cl1),
                               cell_df2 %>% select(accession, cl)
                          )) %>% unique()



# map the cell line field
mapped_sample0 <- mapTextCl(text_df2, cell_df_nodash)                    
mapped_sample1 <- mapTextCl(text_df2 %>% mutate(str=orig_str), 
                            cell_df2 %>% mutate(cl=orig_cl))
mapped_sample2 <- mapTextCl(text_df2 %>% mutate(str=str1), 
                            cell_df_nodash)
mapped_sample3 <- mapTextCl(text_df2 %>% mutate(str=str2), 
                            cell_df_nodash)
mapped_samples <- do.call(rbind, list(mapped_sample0, mapped_sample1, mapped_sample2, 
                                      mapped_sample3)) %>%
                            select(accession, orig_str) %>%
                            unique() # 4359 in human compendia, 1526 in human rnaseq, 511 in mouse compendia, 246 mouse rnaseq
mapped_samples2 <- mapped_samples %>% 
  group_by(orig_str) %>% 
  summarise(accession=paste(accession, collapse=";"))
no_map <- text_df2 %>% select(orig_str) %>% unique() %>% anti_join(mapped_samples)  # 1273, 702 in human, 1473 in mouse
no_map2 <- no_map %>% filter(nchar(orig_str)>3) # 1216, 1402 in mouse
#mapped_samples %>% filter(nchar(orig_str)==3) 
text_mapped <- text_df %>% select(gsm, orig_str) %>% inner_join(mapped_samples2, by=c("orig_str"))

# map by title and part
mapped_sample1_t <- mapTextCl(text_df_title2 %>% mutate(str=orig_str), 
                            cell_df2 %>% mutate(cl=orig_cl), three_l=FALSE)
mapped_sample2_t <- mapTextCl(text_df_title2 %>% mutate(str=str1), 
                            cell_df_nodash, three_l=FALSE)

mapped_samples_t <- do.call(rbind, list( mapped_sample1_t, mapped_sample2_t)) %>%
  select(accession, orig_str) %>%
  unique() %>% 
  group_by(orig_str) %>% 
  summarise(accession=paste(accession, collapse=";"))
t_mapped <- text_df_title %>% inner_join(mapped_samples_t, by=c("orig_str")) %>% unique()

t_mapped_new <- t_mapped %>% filter(!gsm %in% text_mapped$gsm)
t_mapped_new %>% select(orig_str, accession) %>% unique() %>% group_by(accession) %>% 
  mutate(n=n()) %>% sample_n(1) %>% arrange(desc(n)) 
# write this out

length(unique(text_mapped$gsm))==nrow(text_mapped)
length(unique(text_df$gsm))==nrow(text_df)
nrow(text_mapped)/nrow(text_df) # 87.3% mapping!!! compendia, 77.9% rnaseq, 11% mouse compendia, 3% rnaseq
text_mapped %>% write_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl.csv", prefix, data_type))
t_mapped_new %>% write_csv(sprintf("data/02_labeled_data/%s_%s_sample_cl_part2.csv", prefix, data_type))

# add the title mapped ones

#### STUDY MAPPING ####
mapped_study0 <- mapTextCl(study_text_df, cell_df_nodash, three_l=FALSE) 
mapped_study1 <- mapTextCl(study_text_df %>% mutate(str=orig_str), 
                            cell_df2 %>% mutate(cl=orig_cl), three_l=FALSE)
mapped_study2 <- mapTextCl(study_text_df %>% mutate(str=str1), 
                            cell_df_nodash, three_l=FALSE)
mapped_studies <- do.call(rbind, list(mapped_study1, mapped_study2, mapped_study0)) %>%
  select(accession, orig_str) %>%
  unique() #90
mapped_studies2 <- mapped_studies %>% 
  group_by(orig_str) %>% 
  summarise(accession=paste(accession, collapse=";")) # 79
study_to_cl <- study_text_df %>% 
  select(gse, orig_str) %>% 
  inner_join(mapped_studies2, by="orig_str")
study_to_cl %>% write_csv(sprintf("data/02_labeled_data/%s_%s_study_cl.csv", prefix, data_type))

