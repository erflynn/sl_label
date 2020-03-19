require("tidyverse")

list.studies <- list.files(sprintf("data/rnaseq/%s/", prefix), pattern="*RP*")

print(length(list.studies))
list.studies <- list.studies[list.studies != "README.md"]
small.list <- extractChunk(list.studies, idx, SIZE.CHUNK=500)

meta_res <- lapply(list.studies, function(study_id){
  my.f <- sprintf("data/rnaseq/%s/%s/metadata_%s.tsv", prefix, study_id, study_id)
  if (!file.exists(my.f)){
    return(NA)
  }
  meta_ind <- read_tsv(my.f)
  return(meta_ind)
})
meta_df <- do.call(rbind, meta_res[!is.na(meta_res)])
meta_df %>% write_csv("data/human_metadata_tsv_combined.tsv")

# load all the metadata and compare to the main metadata file

metadata <- read.csv("data/01_metadata/human_rnaseq_sample_metadata.csv")
mapping <- read.csv("data/01_metadata/human_rnaseq_exp_to_sample.csv")

setdiff(metadata$acc, meta_df$refinebio_accession_code)
# ok so these are the same -- but the mapping file does overdo things, but that's ok.

# ------------------------- #

load(sprintf("data/%s_counts_rna.RData", prefix))
summary(sapply(res, function(x) x$nrow))
#mode <-- 189440
#2nd most common <-- 185401
#num_rows[num_rows <  183233]
#[1] 183233 145870 107548
# <--- this doesn't show within study NAs :( ---> #

samples <- unique(unlist(sapply(res, function(x) x$samples)))
samples2 <- samples[samples!="gene_name"]
sample_metadata <- read.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", prefix) )
sample_metadata2 <- sample_metadata %>% mutate(present=(acc %in% samples2)) 
table(sample_metadata2$present) # 122864 out of 229789
# mouse:  125652 out of 359142
mapping <- read_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample.csv", prefix))
mapping2 <- mapping %>% mutate(present=(sample_acc %in% samples2))
mapping2 %>% write_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample2.csv", prefix))
sample_metadata2 %>% write.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata2.csv", prefix))

mapping2 %>% filter(present) %>% select(study_acc) %>% unique() %>% count() # 6290 human, 6499 mouse
mapping2 %>% select(study_acc) %>% unique() %>% count() # 6457 human, 6651 mouse

