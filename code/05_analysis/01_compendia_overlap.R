# look at the overlap btw the compendia
# - breakdown into microarray vs RNA-seq

require('tidyverse')
options(stringsAsFactors = FALSE)

compendia <- read.csv("data/01_metadata/human_metadata.csv") # 430119
rnaseq <- read.csv("data/01_metadata/human_rnaseq_sample_metadata.csv") # 229789

overlap_ids <- intersect(compendia$acc, rnaseq$acc) # 99611 overlapping IDs
compendia_only <- setdiff(compendia$acc, rnaseq$acc) # 330508
rnaseq_only <- setdiff(rnaseq$acc, compendia$acc) # 130178
unlist(unique(str_extract_all(rnaseq_only, "[A-Z]+")))
unlist(unique(str_extract_all(overlap_ids, "[A-Z]+")))
sort(unique(unlist(str_extract_all(compendia_only, "[A-Z]+"))))

# *the DRR, SRR are all RNAseq*


# ----- SANITY CHECKS ----- #
# 1. do these have the same sex labels?
compendia_sl <- read_csv("data/02_labeled_data/human_all_sl.csv") %>% 
  mutate(sex=ifelse(pred > 0.5, 1, 0)) # 430,119
rnaseq_sl <- read_csv("data/02_labeled_data/human_rnaseq_sl.csv") %>% 
  mutate(sex=ifelse(pred > 0.5, 1, 0)) # 116,716
overlap_sl <- compendia_sl %>% inner_join(rnaseq_sl, by=c("id")) # 69,852
length(setdiff(rnaseq_sl$id, compendia_sl$id)) # 46,864
# are any of the missing RNA-seq ids in compendia?
compendia_sl %>% filter(id %in% overlap_ids) %>% nrow() # 99,611 -- YUP
# we can add another 30k to the RNA-seq data...

overlap_sl %>% filter(sex.x != sex.y) %>% nrow() # 4087 do not match
# this is 5.9%
overlap_sl %>% filter(sex.x == sex.y & (pred.x < 0.2 | pred.x > 0.8) & 
                        (pred.y < 0.2 | pred.y > 0.8)) %>% nrow() # 33672
overlap_sl %>% filter(sex.x != sex.y & (pred.x < 0.2 | pred.x > 0.8) & 
                        (pred.y < 0.2 | pred.y > 0.8))  
# 244 have high confidence swaps!! jeezus (but this is 0.7%)

# how do these compare to the metadata?
compendia_sex <- read.csv("data/01_metadata/human_microarray_metadata_sex.csv") %>% unique()
rnaseq_sex <- read.csv("data/01_metadata/human_rnaseq_metadata_sex.csv") %>% unique()
metadata_sex <- compendia_sex %>% inner_join(rnaseq_sex, by=c("acc")) # 99611
metadata_sex %>% filter(mapped_sex.x != mapped_sex.y) # none
metadata_lab <- metadata_sex %>% filter(mapped_sex.x %in% c("male", "female")) # 11573
metadata_lab2 <- metadata_lab %>% select(acc, mapped_sex.x) %>% rename(metadata_sex=mapped_sex.x) %>%
  mutate(metadata_sex=ifelse(metadata_sex=="female", 0, 1))

metadata_sl <- metadata_lab2 %>% inner_join(overlap_sl, by=c("acc"="id")) # 7711
metadata_sl %>% filter(metadata_sex!=sex.x ) %>% nrow() # 770 (10%)
metadata_sl %>% filter(metadata_sex!=sex.y ) %>% nrow() # 714 (9.3%)

# ~20% of these are mislabeled in both
metadata_sl %>% filter(metadata_sex!=sex.x & metadata_sex==sex.y) %>% nrow() # 178
metadata_sl %>% filter(metadata_sex!=sex.y & metadata_sex==sex.x) %>% nrow() # 122

# at high confidence predictions, the majority of the mislabels are in both 
metadata_sl %>% filter(metadata_sex != sex.x & (pred.x < 0.2 | pred.x > 0.8)) %>% nrow() # 428
metadata_sl %>% filter(metadata_sex != sex.y & (pred.y < 0.2 | pred.y > 0.8)) %>% nrow() # 162
metadata_sl %>% filter(metadata_sex != sex.x & (pred.x < 0.2 | pred.x > 0.8) & metadata_sex != sex.y) %>% nrow() # 378
metadata_sl %>% filter(metadata_sex != sex.y & (pred.y < 0.2 | pred.y > 0.8) & metadata_sex != sex.x) %>% nrow() # 159

table(metadata_sl$metadata_sex, metadata_sl$sex.x ) 
# errors seem to be more more actually "male", called "female" (592 vs 178)

table(metadata_sl$metadata_sex, metadata_sl$sex.y ) 
# errors seem to be more more actually "male", called "female" (602 vs 112)

# is this a cell line artifact? or something else? I am not clear..

# for the rnaseq data that is only in the compendia, there are few errors (2.4%)
rnaseq_sl_miss <- metadata_lab2 %>% 
  anti_join(overlap_sl, by=c("acc"="id")) %>% 
  inner_join(compendia_sl, by=c("acc"="id")) # 3862
rnaseq_sl_miss %>% filter(metadata_sex != sex) %>% nrow() # 94 (2.4%)
table(rnaseq_sl_miss$metadata_sex,rnaseq_sl_miss$sex ) 
# <-- I do not understand why there are fewer errors here? and less class imbalance?



# 2. do these have the same cell line labels?
compendia_cl <- read_csv("data/02_labeled_data/human_compendia_sample_cl.csv")
rnaseq_cl <- read_csv("data/02_labeled_data/human_rnaseq_sample_cl.csv")
overlap_cl <- compendia_cl %>% inner_join(rnaseq_cl, by=c("gsm")) # 12,227
overlap_cl %>% filter(accession.x != accession.y) 
# all the same, and the same set of overlapping ids! 
