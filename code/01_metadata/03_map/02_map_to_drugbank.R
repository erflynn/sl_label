require('tidyverse')
require('tidytext')
options(stringsAsFactors = FALSE)

# read in the metadata
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
metadata <- read.csv(sprintf("data/%s_metadata2.csv", prefix))

trt_words <- metadata %>% group_by(trt) %>% count()
trt_words2 <- trt_words %>% filter(! trt %in% c("","-", "--", ".") )
#refine_trt_annot <- 
gsm_data <- metadata %>% filter(! trt %in% c("","-", "--")) %>% select(acc, trt) %>% rename(gsm=acc, str=trt)
metadata %>% filter(! trt %in% c("", "--", "-", ".", "n/a") & 
                      is.na(as.numeric(trt)) & 
                      nchar(trt)>3) %>% nrow()
# -- missing a lot of IL

# ----- drug data ------ #
# create the synonym DF
drug_full_info <- read.delim("../labeling/geo2drug/data/00_db_data/drugbank_parsed.txt")
drug_info_df <- read.delim("../labeling/geo2drug/data/00_db_data/drugbank_vocab_no_nutra.txt")

# only match drug data with more than 3 characters and 1 word 
#  TODO - what about drugs that are phrases?
drug_info_df$numchar <- sapply(drug_info_df$name, nchar)
short_names <- filter(drug_info_df, numchar <= 3)
names_to_match <- filter(drug_info_df, numchar > 3)
drug_info_df$nwords <- sapply(drug_info_df$name, function(x) length(strsplit(x, " ")[[1]]))

# fix tabs and trailing punctuation
text_df <- gsm_data[,c("gsm", "str")]
text_df$str <- sapply(text_df$str, function(x)
{y <- gsub(";\t", " ", x);  ## fix ;\t
z <- gsub(' [[:punct:]]|[[:punct:]] ', ' ', y);  ##remove trailing punctuation
return(z)})

# separate out the gsm data into unigrams
gsm_unigrams <- text_df %>% unnest_tokens(word, str, token=stringr::str_split, pattern = " ") 
gsm_unigrams2 <- filter(gsm_unigrams, str_length(word) > 3 & is.na(as.numeric(word))) # remove short words or numbers 

# map with a join
comb_names <- inner_join(filter(drug_info_df, nwords==1), gsm_unigrams2, by=c("name"="word")) %>% distinct()
length(unique(comb_names$gsm)) 
drugs <- unique(comb_names$name ) 
gsm_map_counts <- gsm_unigrams %>% filter (word %in% drugs) %>% group_by(word) %>% summarise(total=n())  %>% arrange(desc(total))


# try some bigrams
bigrams <- gsm_data[,c("gsm", "str")] %>% unnest_tokens(bigram, str, token = "ngrams",  n = 2)
trigrams <- gsm_data[,c("gsm", "str")] %>% unnest_tokens(trigram, str, token = "ngrams", n = 3)

comb_names_bi <- inner_join(filter(drug_info_df, nwords==2), bigrams, by=c("name"="bigram")) %>% distinct()
comb_names_tri <- inner_join(filter(drug_info_df, nwords==3), trigrams, by=c("name"="trigram")) %>% distinct()

#  problem - discarding aa from the drugbank classification
drug_stopwords <- c("same", "dmso", "water", "sage", "camp", "balance",  "biotin")
# need to keep in ATRA, SAHA

# pull in additional curated stopwords
jake_stopwords <- read.delim("../labeling/geo2drug/data/00_db_data/jake_stopwords.txt", head=FALSE)$V1
comb_names2 <- filter(select(comb_names, c("name", "dbID", "gsm")), ! name  %in% union(drug_stopwords, jake_stopwords))
length(unique(comb_names2$gsm)) # 9735
length(unique(comb_names2$dbID)) # 1104


# join with ATC
comb_names3 <- left_join(comb_names2, drug_full_info[,c("dbID", "ATC")], by="dbID")

# NOW TO DISAMBIGUATE - currently has synonymns as names
drug_data_gsm <- inner_join(drug_full_info[,c("dbID", "name")], select(comb_names3, c("dbID", "ATC", "gsm")))


drug2 <- filter(drug_data_gsm, ! name %in% 
                  c("Glucose", "D-glucose", "Oxygen", "Nitrogen", "L-Glutamine", "Sucrose"))
length(unique(drug2$dbID)) # [1] 1098, 253
length(unique(drug2$gsm)) # [1] 8834, 6917 --- out of 39915 samples [34/100 are drugs -- 12906 are ctls] - 32%
# all but 88 of these have study level annotations

write_csv(drug2, sprintf("data/%s_sample_drug_mapped.csv", prefix))
# what is the overlap with study-level annot?
# -- for now, I only have study level annot for mouse/rat in GEO

gse_drug <- read.delim("../labeling/geo2drug/data/02_labeled_data/drugbank_mapped_gse.txt")
human_map <- read_csv("data/human_exp_to_sample.csv") %>% rename(gse=study_acc, gsm=sample_acc)
# how good are these?

drug3 <- drug2 %>% select(gsm, dbID, name) %>% 
  left_join(human_map) %>% 
  full_join(gse_drug %>% select(dbID, name, gse), by="gse")

drug_gse_only <- drug3 %>% filter(str_detect(gse, "GSE"))  

drug_gse_only %>% 
  group_by(gse) %>% filter(!is.na(dbID.x)) %>% count() %>% arrange(n)
# 432 out of 8922 gse have gsm annot

drug_gse_only %>% 
  group_by(gse) %>% filter(!is.na(dbID.y)) %>% count() %>% nrow()
# 8834 have gse annot out of 8922

# how many are the same
present.both <- drug_gse_only %>% filter(!is.na(dbID.x) & !is.na(dbID.y))
present.both %>% filter(dbID.x == dbID.y) %>% nrow() # 5519
matching <- present.both %>% filter(dbID.x == dbID.y)
length(unique(matching$gse)) # 329

missing_gse <- drug_gse_only %>% 
  filter(!is.na(dbID.x)) %>%
  filter(!gse %in% matching$gse) # 103 gses, 1386 gsms
# // THESE are pretty convincing

# ---> but what about ones with GSE and not GSM??
drug_gse_only %>% 
  group_by(gse, name.y) %>% filter(all(is.na(dbID.x))) %>%
  count()
# // some of these are but some of these actually -ARE- trt


# --- what about the 3char ones? --- #
  # these are a non-issue, but we need to look at DASHES
  
# why do we only get 253 instead of 1098?

# what is the overlap w CREEDS



# ------------------------ STOP here --------------------- #


# add a mouse or human flag
drug2_gsm <- left_join(drug2, gsm_data %>% select(gsm, gpl, organism, study_type))

# FYI these are multimapping!
write.table(drug2_gsm, file="data/02_labeled_data/drugbank_mapped_gsm.txt", row.names=FALSE, sep="\t")

d2 <- drug2_gsm %>% select(gsm, organism, study_type) %>% unique()
d2 %>% group_by(organism, study_type) %>% count()
# human    oligo       3883
# human    seq         1017
# mouse    oligo       2862
# mouse    seq         1072

# write out download lists
drug2_gsm %>% filter(organism=="mouse" & study_type=="oligo") %>% select(gsm) %>% unique() %>%
  write_csv("data/01_sample_lists/mouse_oligo_drug.csv") # 2862
drug2_gsm %>% filter(organism=="human" & study_type=="oligo") %>% select(gsm) %>% unique() %>% 
  write_csv("data/01_sample_lists/human_oligo_drug.csv") # 3883

