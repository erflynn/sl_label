

# map both study + sample
cell_df <- read_csv("data/00_db_data/cell_syn_df.csv")
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
# we should be doing on the level of unique tokens

mapped_sample <- mapTextCl(text_df, cell_df)
no_map <- text_df %>% anti_join(mapped_sample)

cell_name <- no_map %>% filter(str_detect(str, "cell"))
cell_line <- cell_name %>% filter(str_detect(str, "line"))
comb_names_n3.2 <- comb_names_n3 %>% filter(!cl %in% char3_stop)
map_cell_3 <- comb_names_n3.2 %>% filter(str_detect(orig_str, "cell"))

# this is not good... -- this means the same synonym maps to *MULTIPLE* cells
map_cell_3 %>% select(cl, accession) %>% unique() %>% arrange(cl, accession) %>% View()

# ---- which contain ATCCC ? map those ---- #
# atcc, crl-1427
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