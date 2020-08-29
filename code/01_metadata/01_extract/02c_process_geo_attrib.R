# Pull sample/attr data from GEOMetadb
#
require('GEOmetadb')
require('tidyverse')

con = dbConnect(SQLite(), "../labeling/GEOmetadb.sqlite")
all_gses = dbGetQuery(con, "SELECT gse FROM gse;")
sample_list = comb_metadata %>% 
  filter(str_detect(sample_acc, "GSM")) %>% 
  #sample_n(20000) %>% 
  pull(sample_acc)

my_gses = by_study %>% filter( str_detect(study_acc, "GSE")) %>% distinct(study_acc)
setdiff(my_gses$study_acc, all_gses$gse) # yayy - all there
dbListFields(con, "gsm")
gsm_dat = dbGetQuery(con, sprintf("SELECT gsm, source_name_ch1, characteristics_ch1, \
                     treatment_protocol_ch1, description FROM gsm \
                     WHERE gsm IN ('%s');", paste(sample_list, collapse="\',\'")))


geometadb_dat <- gsm_dat  %>% as_tibble()
rm(gsm_dat)
save(geometadb_dat, file="data/gsm_meta_geometadb.RData")

separated=geometadb_dat %>% 
  filter(!is.na(characteristics_ch1)) %>% 
  select(gsm, characteristics_ch1)  %>% 
  #sample_n(10000) %>%
  separate_rows(characteristics_ch1, sep=";\t") 
kv <- separated %>%  
  separate(characteristics_ch1, into=c("key", "value"), sep=":", #) %>%
           extra="merge", fill="left") %>%
  mutate(key=ifelse(is.na(key), "other", str_trim(tolower(key))))

# clean up the values
kv %>% write_csv("data/gsm_key_value.csv")

