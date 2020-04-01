
# 0. Remove synonyms
# 1. MGI list of mouse and rat strains
# 2. Edit distance
# 3. 3 character

options(stringsAsFactors = FALSE)
# map both study + sample
cell_df <- read_csv("data/00_db_data/cell_syn_df.csv") # 244810, 109447 acc
n2 <- cell_df %>% filter(numchar < 3) # 379
num <- cell_df %>% filter(!is.na(as.numeric(cl))) # 1222

cell_df2 <- cell_df %>% filter(numchar >= 3 & is.na(as.numeric(cl)))

prefix <- "human"

sample_metadata <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", prefix))
study_metadata <- read_csv(sprintf("data/01_metadata/%s_experiment_metadata.csv", prefix))

refine_cl_annot <- sample_metadata %>% select(acc, cl_line, part, title)  %>%
  filter(!is.na(cl_line) & cl_line != "" & cl_line != "--")
refine_cl_annot %>% filter(str_detect(cl_line, "primary")) %>% select(cl_line) %>% unique()
refine_cl_annot %>% filter(str_detect(cl_line, "stem")) %>% select(cl_line) %>% unique()
refine_cl_annot %>% filter(str_detect(cl_line, "culture")) %>% select(cl_line) %>% unique()



# clean up the input
text_df <- refine_cl_annot %>% rename(gsm=acc, orig_str=cl_line) %>%
  mutate(orig_str=tolower(orig_str)) %>%
  mutate(str=str_trim(str_squish(str_replace_all(orig_str, "[-|,|\\.|/|#|_|(|)|+|;|\t]", ""))))


text_df2 <- text_df %>%
  filter(nchar(str) >= 3 & is.na(as.numeric(str))) 
text_df3 <- text_df2 %>% select(orig_str, str) %>% unique()
#cell_df <- cell_df2

cell_df <- rbind(cell_lab_name, cell_lab_syn %>% rename(cl=synonyms)) %>%
  mutate(
    nwords=length(strsplit(cl, " ")[[1]]),
    numchar=nchar(cl)) 

mapped_sample1 <- mapTextCl(text_df3 %>% mutate(str=orig_str), cell_df)
mapped_sample2 <- mapTextCl(text_df3, cell_df)
mapped_sample3 <- mapTextCl(text_df3, cell_df %>%
                              mutate(cl=str_replace_all(cl, "[-|,|\\.|/|#|_|(|)|+|;|\t]", "")))

mapped_samples <- do.call(rbind, list(mapped_sample1, mapped_sample2, mapped_sample3)) %>%
                            select(accession, orig_str) %>%
                            unique()

no_map <- text_df3 %>% anti_join(mapped_samples)
no_map %>% filter(nchar(str)==3) %>% nrow() # 226
no_map2 <- no_map %>% filter(nchar(str)>3) 

# 3grams

gsm_unigrams <- text_df3 %>%
  unnest_tokens(word, str, token=stringr::str_split, pattern = " ") %>%
  filter( is.na(as.numeric(word))) # remove numbers
gsm_unigrams_char3 <- filter(gsm_unigrams, 
                             str_length(word)==3)

comb_names_n3 <- inner_join(cell_df %>% 
                              mutate(cl=str_replace_all(cl, "[-|,|\\.|/|#|_|(|)|+|;|\t]", "")) %>% 
                              filter(nwords==1), 
                            gsm_unigrams_char3, by=c("cl"="word")) %>% 
  distinct()
gsm_unigrams_char3 %>% group_by(word) %>% count() %>% arrange(desc(n))
char3_stop <- c("the", "and", "sum", "age", "hey", "for", "wbs", "mel", "ctr")
comb_names_n3.2 <- comb_names_n3 %>% filter(!cl %in% char3_stop)
mapped3 <- comb_names_n3.2 %>% select(accession, orig_str) %>% unique() %>%
  group_by(orig_str) %>% summarise(accession=paste(accession, collapse=";"))
mapped3_new <- mapped3 %>% anti_join(mapped_samples, by="orig_str") # 267 new

#### separate things out
mapped_nbsp <- mapTextCl(text_df3 %>% mutate(str=str_replace_all(str, " ", "")), 
                         cell_df %>%
                           mutate(cl=str_replace_all(cl, "[-|,|\\.|/| |#|_|(|)|+|;|\t]", "")))
mapped_nbsp2 <- mapped_nbsp %>% anti_join(mapped_samples, by="orig_str")


mapped_samples2 <- mapped_samples %>% bind_rows(mapped3_new, mapped_nbsp2) 
no_map3 <- no_map2 %>% anti_join(mapped_samples2, by=c("orig_str"))

no_map3 %>% View()
# does a substring exactly match?

# PC3

# CRL



#### EDIT DISTANCE
my.str <- "mda-mb-231-s"
cl.list <- cell_df$cl 


my_fuzzy_match <- function(my.str, cl.list){
  close.matches <- agrep(my.str, cl.list, max.distance=c("insertions"=0.1, "deletions"=0.1, "substitutions"=0), value=TRUE)
  if (length(close.matches) > 1){
    dists <- stringdist::stringdist(my.str, close.matches, method="jw")
    mapping <- close.matches[which(dists==min(dists))] # UGH gives first minimum
    return(mapping)
  } else {
    return(NA)
  }
}
res <- sapply(no_map3$str, function(my.str) my_fuzzy_match(my.str, cl.list))
length(res[!is.na(res)])
# skip the fuzzy mapping

cell_name <- no_map3 %>% filter(str_detect(str, "cell"))
cell_line <- cell_name %>% filter(str_detect(str, "line"))
p_i_s <- no_map3 %>% filter(str_detect(str, "primary|stem|culture|ipsc")) #.. not really sure what this means
# ---- which contain ATCCC ? map those ---- #
# atcc, crl-1427
atcc <- no_map3 %>% filter(stringr::str_detect(orig_str, "atcc|crl"))


# other todos:
#  0. get ATCCC mapping from cellosaurus
#  1. look at "part" & "title"
#  2. is there a cell line within the sentence? grab it if there is




# --- look at 100 that map and 100 that don't --- #
no_map <- sample_metadata %>% 
  anti_join(comb_names , by=c("acc"="gsm")) %>% 
  filter(!is.na(cl_line)) %>%
  group_by(cl_line) %>% 
  count() %>%
  arrange(desc(n))

